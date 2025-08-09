import pyvectorem as em

mesh = em.load_gmsh("examples/cube_cavity.msh")
bc = em.build_scalar_pec(mesh, 1)
params = em.ScalarHelmholtzParams()
params.k0 = 1.0
asm = em.assemble_scalar_helmholtz(mesh, params, bc)
res = em.solve_linear(asm.A, asm.b, em.SolveOptions())
em.write_touchstone("smoke.s2p", [1.0], [em.SParams2()])
print("ok")
