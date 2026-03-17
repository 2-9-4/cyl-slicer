#!/usr/bin/env python3
"""
Interactive G-code viewer for cylindrical slicer output using Plotly.
Proper 3D depth buffering with WebGL. Layer-by-layer navigation.
"""

import numpy as np
import plotly.graph_objects as go
from collections import defaultdict
import re

def parse_gcode(filename):
    """Parse G-code file and extract layer data with section types."""
    layers = defaultdict(lambda: {"moves": [], "z": 0.0})
    current_layer = 0
    current_z = 0.0
    current_type = "UNKNOWN"
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(";LAYER:"):
                current_layer = int(line.split(":")[1])
            elif line.startswith(";Z:"):
                current_z = float(line.split(":")[1])
                layers[current_layer]["z"] = current_z
            elif line.startswith(";TYPE:"):
                current_type = line.split(":")[1]
            elif line.startswith(";"):
                continue
            else:
                parts = line.split()
                if not parts:
                    continue
                cmd = parts[0]
                if cmd in ("G0", "G1"):
                    params = {}
                    for p in parts[1:]:
                        if p[0] in "XYZEF":
                            try:
                                params[p[0]] = float(p[1:])
                            except ValueError:
                                pass
                    layers[current_layer]["moves"].append({
                        "type": "rapid" if cmd == "G0" else "extrude",
                        "section_type": current_type,
                        "params": params
                    })
    
    return dict(layers)


def cylindrical_to_cartesian(theta_deg, x, radius):
    """Convert cylindrical coords to 3D Cartesian.
    
    Coordinate system:
    - theta_deg: angle in degrees (Y in G-code) - rotation around cylinder
    - x: axial position (X in G-code) - along cylinder length
    - radius: cylinder radius
    
    Returns (x, y, z) where X is axial, YZ form the circular cross-section.
    """
    theta = np.radians(theta_deg % 360)
    y = radius * np.cos(theta)
    z = radius * np.sin(theta)
    return x, y, z  # X stays axial


def extract_paths(moves, radius):
    """Extract continuous paths from G-code moves, preserving section types."""
    paths = {"rapid": [], "extrude": []}
    current_path = []
    last_pos = None
    current_section = "UNKNOWN"

    for move in moves:
        params = move["params"]
        if "X" not in params or "Y" not in params:
            continue

        x = params["X"]
        theta = params["Y"]
        pos = (theta, x)
        move_type = move["type"]
        section_type = move.get("section_type", "UNKNOWN")

        if move_type == "rapid":
            if len(current_path) > 1:
                paths["extrude"].append((current_path, current_section))
            current_path = []
            if last_pos is not None:
                current_path.append(last_pos)
            current_path.append(pos)
            paths["rapid"].append((current_path, section_type))
            current_path = []
            current_section = section_type
        else:
            # Extrude move - track section type changes
            if section_type != current_section and current_path:
                # Save current path with its section type
                paths["extrude"].append((current_path, current_section))
                current_path = []
            if not current_path:
                current_section = section_type
                if last_pos is not None:
                    current_path.append(last_pos)
                else:
                    current_path.append(pos)
            current_path.append(pos)
        
        last_pos = pos

    if len(current_path) > 1:
        paths["extrude"].append((current_path, current_section))

    return paths


class GcodeViewer:
    def __init__(self, gcode_file, cylinder_radius=50):
        self.gcode_file = gcode_file
        self.radius = cylinder_radius
        self.layers = parse_gcode(gcode_file)
        self.layer_keys = sorted(self.layers.keys())
        self.current_idx = 0
        
        self.fig = go.Figure()
        self.setup_figure()
        self.update_display()
        
    def setup_figure(self):
        self.fig.update_layout(
            scene=dict(
                aspectmode='data',
                xaxis_title='X (axial)',
                yaxis_title='Y',
                zaxis_title='Z',
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=0.8)
                )
            ),
            title_text='Layer 0 | Z=0.00mm',
            width=1200,
            height=800,
            updatemenus=[{
                'buttons': [
                    {
                        'label': '◀ Prev',
                        'method': 'update',
                        'args': [{}]
                    },
                    {
                        'label': 'Next ▶',
                        'method': 'update', 
                        'args': [{}]
                    }
                ],
                'direction': 'left',
                'x': 0.5,
                'y': 0.02,
                'type': 'buttons',
                'showactive': False
            }]
        )
        
        # Keyboard callbacks via dash would be ideal but keeping it simple
        print("Controls: Use arrow keys in plot window, or click legend to toggle")
        
    def update_display(self):
        self.fig.data = []
        
        layer_num = self.layer_keys[self.current_idx]
        layer_data = self.layers[layer_num]
        moves = layer_data["moves"]
        z_height = layer_data["z"]
        current_radius = z_height
        
        paths = extract_paths(moves, current_radius)
        
        # Compute bounds
        all_x = [m["params"].get("X", 0) for l in self.layers.values() for m in l["moves"]]
        x_min, x_max = min(all_x, default=0), max(all_x, default=100)
        
        # Draw semi-transparent cylinder surface for depth perception
        # This is drawn at the print radius so paths appear on the surface
        theta_cyl = np.linspace(0, 2*np.pi, 50)
        x_cyl = np.linspace(x_min, x_max, 20)
        theta_cyl, x_cyl = np.meshgrid(theta_cyl, x_cyl)
        y_cyl = current_radius * np.cos(theta_cyl)
        z_cyl = current_radius * np.sin(theta_cyl)
        
        self.fig.add_trace(go.Surface(
            x=x_cyl, y=y_cyl, z=z_cyl,
            colorscale=[[0, 'rgba(200,200,200,0.1)'], [1, 'rgba(200,200,200,0.1)']],
            opacity=0.1,
            showscale=False,
            hoverinfo='skip',
            showlegend=False,
            lighting=dict(ambient=0.5, diffuse=0.5, roughness=0.5),
            lightposition=dict(x=0, y=0, z=0)
        ))
        
        # Draw opaque inner cylinder to block view through center
        # Slightly smaller radius so it sits behind the print paths
        inner_radius = current_radius * 0.92
        theta_inner = np.linspace(0, 2*np.pi, 40)
        x_inner = np.linspace(x_min, x_max, 10)
        theta_inner, x_inner = np.meshgrid(theta_inner, x_inner)
        y_inner = inner_radius * np.cos(theta_inner)
        z_inner = inner_radius * np.sin(theta_inner)
        
        self.fig.add_trace(go.Surface(
            x=x_inner, y=y_inner, z=z_inner,
            colorscale=[[0, 'rgba(255,255,255,1.0)'], [1, 'rgba(255,255,255,1.0)']],
            opacity=1.0,
            showscale=False,
            hoverinfo='skip',
            showlegend=True,
            legendgroup='Occluder',
            name='Occluder (toggle off to check alignment)',
            lighting=dict(ambient=0.9, diffuse=0.2, roughness=0.9),
            lightposition=dict(x=0, y=0, z=0)
        ))
        
        # Color mapping for section types
        section_colors = {
            "EXTERNAL_PERIMETER": "#00AA00",  # Green
            "INTERNAL_PERIMETER": "#0066CC",  # Blue
            "SOLID_INFILL": "#CC6600",        # Orange
            "INFILL_RADIAL_FULL": "#CC00CC",  # Magenta
            "INFILL_RADIAL": "#FF6600",       # Orange-red
            "INFILL_AXIAL": "#6600CC",        # Purple
            "UNKNOWN": "#888888"              # Gray
        }
        
        # Group paths by section type for legend toggling
        paths_by_section = defaultdict(list)
        for path, section_type in paths["extrude"]:
            paths_by_section[section_type].append(path)
        
        # Draw extrusion paths grouped by section type (for legend toggling)
        for section_type, section_paths in paths_by_section.items():
            color = section_colors.get(section_type, "#888888")
            
            # Combine all paths of this type into single trace
            all_pts = []
            for path in section_paths:
                theta = [p[0] for p in path]
                x = [p[1] for p in path]
                pts = np.array([cylindrical_to_cartesian(t, xi, current_radius)
                               for t, xi in zip(theta, x)])
                all_pts.append(pts)
            
            if all_pts:
                # Draw each path separately but with same legend entry
                for i, pts in enumerate(all_pts):
                    self.fig.add_trace(go.Scatter3d(
                        x=pts[:, 0], y=pts[:, 1], z=pts[:, 2],
                        mode='lines', line=dict(color=color, width=4),
                        opacity=1.0, 
                        showlegend=(i==0),  # Only first path shows in legend
                        legendgroup=section_type,
                        name=section_type.replace('_', ' ')
                    ))
        
        # Draw rapid moves (ALL paths, faint gray) - single trace for legend
        rapid_pts = []
        for path, section_type in paths["rapid"]:
            theta = [p[0] for p in path]
            x = [p[1] for p in path]
            pts = np.array([cylindrical_to_cartesian(t, xi, current_radius)
                           for t, xi in zip(theta, x)])
            rapid_pts.append(pts)
        
        if rapid_pts:
            for i, pts in enumerate(rapid_pts):
                self.fig.add_trace(go.Scatter3d(
                    x=pts[:, 0], y=pts[:, 1], z=pts[:, 2],
                    mode='lines', line=dict(color='red', width=2, dash='dash'),
                    opacity=0.15, 
                    showlegend=(i==0),
                    legendgroup='Travel',
                    name='Travel (rapid)'
                ))
        
        # Update title and layout with legend
        self.fig.update_layout(
            title_text=f'Layer {layer_num} | Z={z_height:.2f}mm | {self.current_idx+1}/{len(self.layer_keys)} | Paths: {len(paths["extrude"])}',
            showlegend=True,
            legend=dict(
                x=1.02,
                y=1,
                xanchor='left',
                yanchor='top',
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1,
                font=dict(size=10)
            ),
            margin=dict(l=0, r=200, t=50, b=0)
        )
        
    def next_layer(self):
        if self.current_idx < len(self.layer_keys) - 1:
            self.current_idx += 1
            self.update_display()
            self.fig.show()
    
    def prev_layer(self):
        if self.current_idx > 0:
            self.current_idx -= 1
            self.update_display()
            self.fig.show()
    
    def show(self):
        self.fig.show()
        
    def run_with_controls(self):
        """Run with keyboard controls using keyboard module."""
        try:
            import keyboard
        except ImportError:
            print("Install keyboard module: pip install keyboard")
            self.fig.show()
            return
            
        print("Controls: Left/Right arrows to navigate, Q to quit")
        
        def on_next(event):
            if event.name == 'right':
                self.next_layer()
            elif event.name == 'left':
                self.prev_layer()
            elif event.name == 'q':
                exit()
        
        keyboard.on_press(on_next)
        self.fig.show()
        keyboard.unhook_all()


def main():
    import sys
    gcode_file = sys.argv[1] if len(sys.argv) > 1 else "output.gcode"
    
    print(f"Loading G-code from: {gcode_file}")
    
    viewer = GcodeViewer(gcode_file, cylinder_radius=50)
    viewer.show()
    
    # Simple text-based navigation
    while True:
        cmd = input("\nCommands: [n]ext, [p]rev, [q]uit, layer number: ").strip()
        if cmd.lower() == 'q':
            break
        elif cmd.lower() == 'n' or cmd == '':
            if viewer.current_idx < len(viewer.layer_keys) - 1:
                viewer.current_idx += 1
                viewer.update_display()
                viewer.fig.show()
        elif cmd.lower() == 'p':
            if viewer.current_idx > 0:
                viewer.current_idx -= 1
                viewer.update_display()
                viewer.fig.show()
        elif cmd.isdigit():
            idx = int(cmd)
            if 0 <= idx < len(viewer.layer_keys):
                viewer.current_idx = idx
                viewer.update_display()
                viewer.fig.show()


if __name__ == "__main__":
    main()
