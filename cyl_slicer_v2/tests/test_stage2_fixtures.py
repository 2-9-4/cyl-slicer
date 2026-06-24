from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np

from cyl_slicer_v2.fixtures import FIXTURES, make_fixture
from cyl_slicer_v2.mesh_io import load_mesh, radial_bounds


class FixtureTests(unittest.TestCase):
    def test_generated_fixtures_are_meshes_with_expected_bounds(self):
        for name in FIXTURES:
            with self.subTest(name=name):
                mesh = make_fixture(name)
                self.assertGreater(len(mesh.vertices), 0)
                r_min, r_max = radial_bounds(mesh)
                self.assertGreater(r_max, r_min)
                self.assertTrue(np.isfinite(mesh.bounds).all())

    def test_watertight_expected_fixtures(self):
        for name in ("solid_box", "tube", "thin_ring", "overhang_step"):
            with self.subTest(name=name):
                self.assertTrue(make_fixture(name).is_watertight)

    def test_existing_repo_stls_load_as_smoke_inputs(self):
        root = Path(__file__).resolve().parents[2]
        for rel in ("tube.stl", "cubeandring.stl", "helix.stl"):
            with self.subTest(stl=rel):
                mesh = load_mesh(root / rel)
                self.assertGreater(len(mesh.faces), 0)


if __name__ == "__main__":
    unittest.main()
