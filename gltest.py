from vispy import gloo, app, scene
app.use_app('pyqt6')  # or pyqt6, pyside2, pyside6
canvas = app.Canvas()
canvas.show()
canvas = scene.SceneCanvas(keys='interactive', show=True)
print("GL Version:", gloo.gl.glGetString(gloo.gl.GL_VERSION).decode())
print("GL Renderer:", gloo.gl.glGetString(gloo.gl.GL_RENDERER).decode())