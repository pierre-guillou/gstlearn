import gstlearn as gl
import numpy as np

dbin = gl.Db.createFromOnePoint(np.zeros(3))
dbin.setLocatorByColIdx(3, gl.ELoc.Z)
dbout = gl.Db.createFromOnePoint(np.ones(2))
model = gl.Model(gl.CovContext())
model.addCov(gl.CovAniso.createIsotropic(model.getContext(), gl.ECov.LINEAR, 1))
ksys = gl.KrigingSystem(dbin, dbout, model, gl.NeighUnique())
ksys.updKrigOptEstim(dbout.addColumnsByConstant(1, 0.0), -1, -1)
ksys.setKrigOptCalcul(gl.EKrigOpt.POINT)
ksys.isReady()  # Undefined behavior, gl.NeighUnique() does not exist anymore
