/***************************************************************************/
/*                                                                         */
/*  Typemap SWIG simple pour gstlrn::ColID&&                               */
/*                                                                         */
/***************************************************************************/
%typemap(in) gstlrn::ColID&& (gstlrn::ColID *temp_colid = nullptr)
{
    bool debug = true; // Passer à false pour désactiver le debug

    SEXP obj = $input;
    SEXP target = obj;

    gstlrn::Id index = 0;
    gstlrn::Id version = 0;

    if (debug) Rprintf("\n[DEBUG Typemap ColID] --- Début conversion ---\n");

    // -------------------------------------------------------------------------
    // 1. Gestion des listes R à 2 éléments : list(objet, param)
    // -------------------------------------------------------------------------
    if (TYPEOF(obj) == VECSXP && Rf_length(obj) == 2)
    {
        target = VECTOR_ELT(obj, 0);
        SEXP param = VECTOR_ELT(obj, 1);
        gstlrn::Id val = Rf_isInteger(param) ? INTEGER(param)[0] : static_cast<gstlrn::Id>(REAL(param)[0]);

        // Seul ERole (objet R ou enum R) utilise le 2ème argument comme 'index'
        if (Rf_inherits(target, "ERole")) {
            index = val;
            if (debug) Rprintf("[DEBUG Typemap ColID] Liste détectée : ERole (index=%d)\n", (int)index);
        } else {
            version = val;
            if (debug) Rprintf("[DEBUG Typemap ColID] Liste détectée : Autre objet (version=%d)\n", (int)version);
        }
    }

    // Extraction du pointeur externe R (slots SWIG S4/Env le cas échéant)
    SEXP ptr_obj = target;
    if (TYPEOF(target) == S4SXP || Rf_isObject(target) || TYPEOF(target) == ENVSXP) {
        SEXP ref = Rf_getAttrib(target, Rf_install("ref"));
        if (ref != R_NilValue) ptr_obj = ref;
    }

    // -------------------------------------------------------------------------
    // 2. Tests successifs des syntaxes de ColID
    // -------------------------------------------------------------------------

    // --- Syntaxe A : Nom de colonne (String) ---
    if (Rf_isString(target) && Rf_length(target) > 0)
    {
        const char* name = CHAR(STRING_ELT(target, 0));
        if (debug) Rprintf("[DEBUG Typemap ColID] Succès -> Par nom: '%s' (version=%d)\n", name, (int)version);
        temp_colid = new gstlrn::ColID(std::string(name), version);
    }

    // --- Syntaxe B : ColID déjà instancié ---
    else if (Rf_inherits(target, "ColID"))
    {
        gstlrn::ColID *col_ptr = nullptr;
        SWIG_ConvertPtr(ptr_obj, (void**)&col_ptr, SWIGTYPE_p_gstlrn__ColID, 0);
        if (col_ptr) {
            if (debug) Rprintf("[DEBUG Typemap ColID] Succès -> Copie de ColID\n");
            temp_colid = new gstlrn::ColID(*col_ptr);
        }
    }

    // --- Syntaxe C : RoleID ---
    else if (Rf_inherits(target, "RoleID"))
    {
        gstlrn::RoleID *roleid_ptr = nullptr;
        SWIG_ConvertPtr(ptr_obj, (void**)&roleid_ptr, SWIGTYPE_p_gstlrn__RoleID, 0);
        if (roleid_ptr) {
            if (debug) Rprintf("[DEBUG Typemap ColID] Succès -> Par RoleID (version=%d)\n", (int)version);
            temp_colid = new gstlrn::ColID(*roleid_ptr, version);
        }
    }

    // --- Syntaxe D : ERole (Objet R ou Enum) ---
    else if (Rf_inherits(target, "ERole"))
    {
        if (TYPEOF(ptr_obj) == EXTPTRSXP) {
            gstlrn::ERole *erole_ptr = nullptr;
            SWIG_ConvertPtr(ptr_obj, (void**)&erole_ptr, SWIGTYPE_p_gstlrn__ERole, 0);
            if (erole_ptr) {
                if (debug) Rprintf("[DEBUG Typemap ColID] Succès -> Par ERole classe (index=%d, version=%d)\n", (int)index, (int)version);
                temp_colid = new gstlrn::ColID(*erole_ptr, index, version);
            }
        } else {
            int enum_val = Rf_isInteger(target) ? INTEGER(target)[0] : static_cast<int>(REAL(target)[0]);
            if (debug) Rprintf("[DEBUG Typemap ColID] Succès -> Par ERole enum=%d (index=%d, version=%d)\n", enum_val, (int)index, (int)version);
            temp_colid = new gstlrn::ColID(static_cast<gstlrn::ERole>(enum_val), index, version);
        }
    }

    // --- Syntaxe E : Numéro d'index brut de colonne (Integer / Real) ---
    else if (Rf_isInteger(target) || Rf_isReal(target))
    {
        gstlrn::Id icol = Rf_isInteger(target) ? static_cast<gstlrn::Id>(INTEGER(target)[0])
                                                : static_cast<gstlrn::Id>(REAL(target)[0]);
        if (debug) Rprintf("[DEBUG Typemap ColID] Succès -> Par index de colonne: %d (version=%d)\n", (int)icol, (int)version);
        temp_colid = new gstlrn::ColID(icol, version);
    }

    // --- Échec ---
    else
    {
        if (debug) Rprintf("[DEBUG Typemap ColID] ÉCHEC -> Type R non reconnu !\n");
        Rf_error("Impossible de convertir l'objet R en ColID");
    }

    $1 = temp_colid;
}

%typemap(freearg) gstlrn::ColID&&
{
    if ($1) {
        delete $1;
    }
}
