
the_ClinVar_Table = 
    ClinVar_Table_CM_SNV[ClinVar_Table_CM_SNV$Name %in% 
                             c('NM_000256.3(MYBPC3):c.2827C>T (p.Arg943Ter)',
                               'NM_001103.4(ACTN2):c.355G>A (p.Ala119Thr)',
                               'NM_170707.4(LMNA):c.673C>T (p.Arg225Ter)'),]
props=theprops
the_editor='Sp_Cas9_Abe8'
    
testTable_out = get_statistics_for_editor(the_editor = the_editor, props = theprops, the_ClinVar_Table = the_ClinVar_Table)
    # testTable_out = testTable

View(testTable_out[testTable_out$Name %in% c('NM_000256.3(MYBPC3):c.2827C>T (p.Arg943Ter)',
                                   'NM_001103.4(ACTN2):c.355G>A (p.Ala119Thr)',
                                   'NM_170707.4(LMNA):c.673C>T (p.Arg225Ter)'),])
