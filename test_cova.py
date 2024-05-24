import gstlearn as gl
import numpy as np


class DualKrigingFunctor:
    def __init__(self, samples, model):
        self._samples = [gl.SpacePoint(p) for p in samples[:, :2]]
        model = gl.Model(model)  # Copy is required in my application
        neigh = gl.NeighUnique()  # Fixes Example 1 seg fault
        self._cov_list = model.getCovAnisoList()
        # Build the Dbs
        dbin = gl.Db.createFromSamples(
            nech=len(samples),
            order=gl.ELoadBy.SAMPLE,
            tab=samples.ravel(),
            names=["x1", "x2", "z1"],
            locatorNames=["x1", "x2", "z1"],
            flag_add_rank=0,
        )
        dbout = gl.Db.createFromOnePoint(np.zeros(2), flag_add_rank=0)
        iptr_est = dbout.addColumnsByConstant(1, 0.0)
        # Build the kriging system
        ksys = gl.KrigingSystem(dbin, dbout, model, neigh)
        ksys.updKrigOptEstim(iptr_est, -1, -1)
        ksys.setKrigOptCalcul(gl.EKrigOpt.POINT)
        ksys.isReady()
        # Solve the kriging system and store the result
        ksys.estimate(0)
        self._zam = np.asarray(ksys.getZam())

    def __call__(self, point):
        spoint = gl.SpacePoint(point)
        c0 = np.array(
            [self._cov_list.eval(s, spoint, 0, 0) for s in self._samples],
            dtype=np.float64,
        )  # Seg fault when calling eval()
        return np.inner(c0, self._zam)


samples = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1]], dtype=np.float64)
model = gl.Model(gl.CovContext())
model.addCov(gl.CovAniso.createIsotropic(model.getContext(), gl.ECov.LINEAR, 1))
functor = DualKrigingFunctor(samples, model)
functor(np.zeros(2))  # Undefined behavior
