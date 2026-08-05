#
# This test is meant to check the manipulation of the new DataBase under Python
#

# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))

######################################################
# We intantiate the DbData object and load information

data = DbData()
invisible(data$printContents("Checking that the DbData is empty"))

invisible(mestitle(0, "Adding Columns of various types"))
writeLines("- Column 0 of type Double with role X, filled with [1.0, 2.0, 3.0]")
invisible(data$addColumnD("hello", VectorDouble(c(1.0, 2.0, 3.0)), RoleID(ERole_X())))

writeLines("- Column 1 of type Int, filled with [5, 6, 7]")
invisible(data$addColumnI("world", VectorInt(c(5, 6, 7)), RoleID(ERole_Z())))

writeLines("- Column 2 of type Bool, filled with [1, 0, 1]")
invisible(data$addColumnB("foobool", VectorBool(c(1, 0, 1))))

writeLines("- Column 3 of type String, filled with ['foo', 'bar', 'baz']")
# Problem: constructing directly the VectorString from a list of strings is not working, so we use the factory
# a = VectorString()
# invisible(a$push_back("foo"))
# invisible(a$push_back("bar"))
# invisible(a$push_back("baz"))
# invisible(data$addColumnS("foobar", a, RoleID(ERole_Z())))

# writeLines("- Column 4 of type Double (5 versions) constantly filled with 3")
# data$addColumnEmptyD("Bonjour", valinit = 3.0, nversion = 5, roleID = RoleID(ERole_F()))

invisible(data$printContents("\nDbData is not empty anymore"))

# ####################################################################
# # In this part, we check the different ways to construct the ColID #

# invisible(mestitle(0, "Testing Various ways to specify a column (e.g. for retreiving its name):"))
# writeLines(paste("- by Column Name:", data$getName("hello")))
# writeLines(paste("- by Column Name and Version:", data$getName(list("Bonjour", 1))))

# writeLines(paste("- by Column Index:", data$getName(4)))
# writeLines(paste("- by Column Index and Version:", data$getName(list(4, 1))))

# colid = ColID(1)
# writeLines(paste("- by ColID:", data$getName(colid)))

# roleid = RoleID(ERole_F())
# writeLines(paste("- by RoleID and Version:", data$getName(list(roleid, 0))))
# writeLines(paste("- by RoleID:", data$getName(roleid)))

# # TODO: probme in next lines
# # role = ERole_F()
# # writeLines(paste("- by Role:", data$getName(role)))
# # writeLines(paste("- by Role and Index:", data$getName(list(role,0))))
# # # Next line is the Unique way to define a ColID based on a Role and specifying the Index and the Version
# # writeLines(paste("- by Role and Index then Version:", data$getName(list(RoleID(role, 0), 2))))

# #####################################################################
# # In this part, we check the different ways to enquiry the DataBase #

# invisible(mestitle(0, "Checking the presence of the columns in the DbData:"))
# writeLines(paste("- Column called 'hello': ", data$hasColumn("hello")))
# writeLines(paste("- Column called 'world': ", data$hasColumn("world")))
# writeLines(paste("- Column called 'foobar': ", data$hasColumn("foobar")))

# invisible(mestitle(0, "Checking that we can retrieve the whole contents of any column:"))

# writeLines(paste0(" - Column 1: (", paste(data$getColumnI(1), collapse = ", "), ")"))
# writeLines(paste0(" - Column 2: (", paste(data$getColumnB(2), collapse = ", "), ")"))
# writeLines(paste0(" - Column 0: (", paste(data$getColumnD(0), collapse = ", "), ")"))
# writeLines(paste0(" - Column 3: (", paste(data$getColumnS(3), collapse = ", "), ")"))

# isample = 2
# invisible(mestitle(0, "Checking that we can manipulate one Target Sample from any column:"))
# writeLines(paste("Retrieving a Target Sample (", isample, ") from any column"))
# writeLines(paste("- From column 0 (specified by its roleID): ", data$getValueD(roleid, isample)))
# writeLines(paste("- From column 1 (specified by its name): ", data$getValueI("world", isample)))
# # TODO: problem in next line
# # writeLines(paste("- From column 2 (specified by its index): ", data$getValueS(2, isample)))

# invisible(mestitle(1, "Modifying the value at Target: Column 0: 4.0, Column 1: 8, Column 2: 'foobar'"))
# data$setValueD(0, isample, 4)
# data$setValueI(1, isample, 8)
# data$setValueS(2, isample, "foobar")

# invisible(mestitle(1, "Checking the new value of the Target Sample"))
# writeLines(paste("Contents of element #", isample, " for Column 0: ", data$getValueD(0, isample)))
# writeLines(paste("Contents of element #", isample, " for Column 1: ", data$getValueI(1, isample)))
# writeLines(paste("Contents of element #", isample, " for Column 2: ", data$getValueS(2, isample)))

# #########################################################
# # In this part, we check the volontary misuse of DbData #

# invisible(mestitle(0, "Misuses of the DbData"))
# data$printContents("- Initial situation")

# writeLines(paste("\nAdding a column with a name that already exists (world) but different type (Double)"))
# data$addColumnD("world", VectorDouble(c(10.0, 11.0, 12.0)))

# writeLines(paste("\nAdding a Column with an already existing Role (X) but non consecutive index (10)"))
# data$addColumnD("hello", VectorDouble(c(101.0, 102.0, 103.0)), RoleID(ERole_X(), 10))

# writeLines(paste("\nAdding a Column with an already existing Role (X) and existing Index (0)"))
# data$addColumnD("hello", VectorDouble(c(101.0, 102.0, 103.0)), RoleID(ERole_X(), 0))

# ##############################
# # Deleting Columns of DbData #

# invisible(mestitle(0, "Deleting a Column (world)"))
# data$printContents("- Initial situation")
# data$deleteColumn("world.1")
# data$printContents("- Final situation")

# ##########################################################
# # Playing with multiple versions in a Column of a DbData #

# invisible(mestitle(0, "Testing the MultiVersion feature of the DbData"))
# data$addColumnEmptyD("Bonjour", nsamples=3, nversion=3, valinit=0.0, forbidNA=TRUE)
# data$printContents("After adding the 'Bonjour' column (with 3 versions)")

# ##################################
# # Additional inquiries on DbData #

# invisible(mestitle(0, "Various inquiries on the DbData"))
# writeLines(paste("- Number of columns: ", data$getNCols()))
# writeLines(paste("- Number of samples: ", data$getNSamples()))
# writeLines(paste("- Number of Versions in column 0: ", data$getNVersions(0)))

# ############################
# # Testing errors on DbData #

# invisible(mestitle(0, "Erroneous operations on the DbData"))
# data$printContents("Initial situation")

# writeLines(paste("\nTrying to delete a non existing column (world.22)"))
# data$deleteColumn("world.22")

# writeLines(paste("\nTrying to add a Column with a different number of samples"))
# data$addColumnD("hello", VectorDouble(c(201.0, 202.0, 203.0, 204.0))))

# writeLines(paste("\nTrying to set a value in a non existing column (world.22)"))
# data$setValueD("world.22", 0, 52.0)

# writeLines(paste("\nTrying to modify the contents of the Column 'hello' (double)"))
# data$setValueS("hello", 0, "Invalid String")

# writeLines(paste("\nTrying to use a wrong version (5) in Column 'Bonjour' (3 versions)"))
# data$getValueD(colID("Bonjour", 5), 0)

# writeLines(paste("\nTrying to set a value to NA in Column 'Bonjour' (forbidNA = TRUE)"))
# data$setValueD("Bonjour", 0, TEST)
