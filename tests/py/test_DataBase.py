import gstlearn as gl
import numpy as np
import pandas as pd

#
# This test is meant to check the manipulation of the new DataBase under Python
#

######################################################
# We intantiate the DbData object and load information

data = gl.DbData()
data.printContents("Checking that the DbData is empty")

gl.mestitle(0, "Adding Columns of various types")
# >>> We should be able to remove converting into gl.VectorDouble() through a dedicated template function

print("- Column 0 of type Double with role X, filled with [1.0, 2.0, 3.0]")
data.addColumnD("hello", gl.VectorDouble([1.0, 2.0, 3.0]), gl.RoleID(gl.ERole.X))

print("- Column 1 of type Int, filled with [5, 6, 7]")
data.addColumnI("world", gl.VectorInt([5, 6, 7]), gl.RoleID(gl.ERole.Z))

print("- Column 2 of type String, filled with ['foo', 'bar', 'baz']")
data.addColumnS("foobar", gl.VectorString(["foo", "bar", "baz"]), gl.RoleID(gl.ERole.Z))

print("- Column 3 of type Bool, filled with [True, False, True]")
data.addColumnB("foobool", gl.VectorBool([True, False, True]))

print("- Column 4 of type Double, with 5 versions and role F, filled with 3.0")
data.addColumnEmptyD("Bonjour", 0, 5, gl.RoleID(gl.ERole.F), 3.0)
data.printContents("\nInitial Data Base contents")

####################################################################
# In this part, we check the different ways to construct the ColID #

gl.mestitle(0, "Various ways to specify a column (e.g. for retreiving its name):")
print("- by Name :", data.getName("hello"))
print("- by Name and Version :", data.getName(("Bonjour", 1)))
print("- by Index :", data.getName(4))
print("- by Index and Version :", data.getName((4, 1)))
colid = gl.ColID(1)
print("- by ColID :", data.getName(colid))
roleid = gl.RoleID(gl.ERole.F)
print("- by RoleID :", data.getName(roleid))
print("- by RoleID and Version :", data.getName((roleid, 0)))
print("- by Role :", data.getName(gl.ERole.Z))
print("- by Role and Index :", data.getName((gl.ERole.Z, 0)))
# Next line is the Unique way to define a ColID based on a Role and specifying the Index and the Version
print("- by Role and Index and Version :", data.getName((gl.RoleID(gl.ERole.Z, 0), 2)))
#####################################################################
# In this part, we check the different ways to enquiry the DataBase #

gl.mestitle(0, "Checking the presence of the columns in the DbData:")
print("- Column called 'hello': ", data.hasColumn("hello"))
print("- Column called 'world': ", data.hasColumn("world"))
print("- Column called 'foobar': ", data.hasColumn("foobar"))

gl.mestitle(0, "Checking that we can retrieve the whole contents of any column:")
print(" - Column 0:", data.getColumnD(0))
print(" - Column 1:", data.getColumnI(1))
print(" - Column 2:", data.getColumnS(2))

isample = 2
gl.mestitle(0, "Checking that we can manipulate one Target Sample from any column:")
print("Retrieving a Target Sample (", isample, ") from any column")
print("- From column 0 (specified by its role): ", data.getValueD(gl.ERole.X, isample))
print("- From column 1 (specified by its name): ", data.getValueI("world", isample))
print("- From column 2 (specified by its index): ", data.getValueS(2, isample))

gl.mestitle(
    1, "Modifying the value at Target: Column 0: 4.0, Column 1: 8, Column 2: 'foobar'"
)
data.setValueD(0, isample, 4)
data.setValueI(1, isample, 8)
data.setValueS(2, isample, "foobar")

gl.mestitle(1, "Checking the new value of the Target Sample")
print("Contents of element #", isample, " for Column 0: ", data.getValueD(0, isample))
print("Contents of element #", isample, " for Column 1: ", data.getValueI(1, isample))
print("Contents of element #", isample, " for Column 2: ", data.getValueS(2, isample))

#########################################################
# In this part, we check the volontary misuse of DbData #

gl.mestitle(0, "Misuses of the DbData")
data.printContents("- Initial situation")

print(
    "\nAdding a column with a name that already exists (world) but different type (Double)"
)
data.addColumnD("world", gl.VectorDouble([10.0, 11.0, 12.0]))

print(
    "\nAdding a Column with an already existing Role (X) but non consecutive index (10)"
)
data.addColumnD(
    "hello", gl.VectorDouble([101.0, 102.0, 103.0]), gl.RoleID(gl.ERole.X, 10)
)

print("\nAdding a Column with an already existing Role (X) and existing Index (0)")
data.addColumnD(
    "hello", gl.VectorDouble([101.0, 102.0, 103.0]), gl.RoleID(gl.ERole.X, 0)
)

##############################
# Removing Columns of DbData #

gl.mestitle(0, "Deleting a Column (world)")
data.printContents("- Initial situation")
data.deleteColumn("world.1")
data.printContents("- Final situation")

##########################################################
# Playing with multiple versions in a Column of a DbData #

gl.mestitle(0, "Testing the MultiVersion feature of the DbData")
data.addColumnEmptyD("Bonjour", nsamples=3, nversion=3, valinit=0.0, forbidNA=True)
data.printContents("After adding the 'Bonjour' column (with 3 versions)")

##################################
# Additional inquiries on DbData #

gl.mestitle(0, "Various inquiries on the DbData")
print("- Number of columns: ", data.getNCols())
print("- Number of samples: ", data.getNSamples())
print("- Number of Versions in column 0: ", data.getNVersions(0))

############################
# Testing errors on DbData #

gl.mestitle(0, "Erroneous operations on the DbData")
data.printContents("Initial situation")

print("\nTrying to delete a non existing column (world.22)")
data.deleteColumn("world.22")

print("\nTrying to add a Column with a different number of samples")
data.addColumnD("hello", gl.VectorDouble([201.0, 202.0, 203.0, 204.0]))

print("\nTrying to set a value in a non existing column (world.22)")
data.setValueD("world.22", 0, 52.0)

print("\nTrying to modify the contents of the Column 'hello' (double)")
data.setValueS("hello", 0, "Invalid String")

print("\nTrying to use a wrong version (5) in Column 'Bonjour' (3 versions)")
colid = gl.ColID("Bonjour", 5)
data.getValueD(("Bonjour", 5), 0)

print("\nTrying to set a value to NA in Column 'Bonjour' (forbidNA = True)")
data.setValueD("Bonjour", 0, gl.TEST)
