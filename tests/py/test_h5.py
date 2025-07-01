import gstlearn as gl

dims = [10, 10, 10]
db = gl.DbGrid.create(nx=dims)
db.dumpToNF("dbgrid.h5", gl.EFormatNF.H5)

db2 = gl.DbGrid.createFromNF("dbgrid.h5")
