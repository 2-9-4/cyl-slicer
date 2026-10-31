import unittest

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.layers import layer_slabs, slice_mesh
from cyl_slicer_v2.toolpaths import gap_support_slices
from cyl_slicer_v2.models import LayerSlice
from shapely.geometry import GeometryCollection, box
import trimesh


def mesh():
    result = trimesh.creation.box(extents=(20, 20, 10))
    result.apply_translation((30, 0, 0))
    return result


class FirstLayerHeightTests(unittest.TestCase):
    def test_cli_profiles_and_validation(self):
        for profile in ('quality', 'legacy-cli'):
            args = build_parser().parse_args(['--profile', profile])
            self.assertIsNone(_profile_value(args, 'first_layer_height'))
            args = build_parser().parse_args(['--profile', profile, '--first-layer-height', '.24'])
            self.assertEqual(_profile_value(args, 'first_layer_height'), .24)
        for height in (0, -1, float('nan'), .5):
            with self.assertRaises(ValueError):
                SlicerConfig(first_layer_height=height, nozzle_diameter=.4).validate()

    def test_first_slab_and_following_slabs(self):
        model = mesh()
        for mode in ('fixed', 'complete-terminal'):
            config = SlicerConfig(bed_radius=30, min_radius=30, max_radius=30.52,
                                  first_layer_height=.2, radial_layer_mode=mode)
            slabs = layer_slabs(model, config)
            for (radius, thickness), (expected_radius, expected_thickness) in zip(
                    slabs, [(30.1, .2), (30.28, .16), (30.44, .16)]):
                self.assertAlmostEqual(radius, expected_radius)
                self.assertAlmostEqual(thickness, expected_thickness)
            self.assertEqual(len(slabs), 3)
            layers = slice_mesh(model, config, max_layers=1)
            self.assertAlmostEqual(layers[0].deposited_thickness, .2)
            self.assertAlmostEqual(layers[0].radius + layers[0].deposited_thickness * .5, 30.2)
            first = config.for_deposited_thickness(layers[0].deposited_thickness)
            self.assertAlmostEqual(first.extrusion_per_mm / config.extrusion_per_mm,
                                   first.bead_area / config.bead_area)

    def test_terminal_remainder_does_not_change_first_height(self):
        config = SlicerConfig(bed_radius=30, min_radius=30, max_radius=30.39,
                              first_layer_height=.2, radial_layer_mode='complete-terminal')
        slabs = layer_slabs(mesh(), config)
        self.assertAlmostEqual(slabs[0][1], .2)
        self.assertAlmostEqual(slabs[1][1], .095)
        self.assertAlmostEqual(slabs[2][1], .095)

    def test_terminal_rejects_unprintably_small_remainder(self):
        config = SlicerConfig(bed_radius=30, min_radius=30, max_radius=30.24,
                              first_layer_height=.2, radial_layer_mode='complete-terminal')
        with self.assertRaisesRegex(ValueError, 'Remaining radial range'):
            layer_slabs(mesh(), config)

    def test_foundation_uses_first_height(self):
        empty = GeometryCollection()
        shape = box(.1, 0, .9, 8)
        layer = LayerSlice(0, 30.45, shape, empty, empty, empty, empty,
                           deposited_thickness=.2)
        for mode in ('fixed', 'complete-terminal'):
            config = SlicerConfig(bed_radius=30, first_layer_height=.2,
                                  support_gap_layers=0, radial_layer_mode=mode)
            support = gap_support_slices([layer], config)
            self.assertAlmostEqual(support[0].deposited_thickness, .2)
            self.assertAlmostEqual(support[0].radius, 30.1)
            self.assertLessEqual(support[-1].radius + support[-1].deposited_thickness / 2,
                                 30.35 + 1e-9)
