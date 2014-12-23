-- Lua script.
p=tetview:new()
p:load_mesh("I:/Programs/VegaFEM-v2.1/myProject/tetGen/tshape.1.ele")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
