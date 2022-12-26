import sys
import hail as hl
hl.init()
################ Usage ##################
# 
# set following arguments.
# path : input_file path
# lcr_bed : bed file require high memory, use lcr bed files splited by chromosome.
# MAF : 0.001, 0.0001, 0.00001
#
# <usage>
# enhancer_promoter_burden_preprocess_Nov182022_hj.py [input_file] [lcr_bed] [MAF]
#
# GQ : 20(dominant), 10(recessive)
# DP : 8(dominant), 4(recessive)
# PLDIFFDP : 8(dominant), NA(recessive)

if len(sys.argv) != 4:
    print("please, check arguments")
else : 
    input_file = sys.argv[1]
    lcr_bed = sys.argv[2]
    MAF = sys.argv[3]
    

class AnnovarWithHail: 
      
    # Called by annovar_vcf_preprocess
    def str_str_preprocess(self, mt, origFieldName): # str -> str
        # hl.is_missing does not work in these columns.
        # because array<str> deal with NA as "NA", the function can not recognize NA.
        obj = hl.if_else(mt.info[origFieldName][0] == 'NA',
                         hl.missing('str'),
                         hl.str(mt.info[origFieldName][0]))
        return obj
    
    def float_float_preprocess(self, mt, origFieldName): # int or float -> float
        obj = hl.if_else(hl.is_missing(mt.info[origFieldName][0]),
                         hl.missing('float64'),
                         hl.float64(mt.info[origFieldName][0]))
        return obj
    
    # Called by annovar_vcf_preprocess
    def str_float_preprocess(self, mt, origFieldName) : # str -> float
        # hl.is_missing does not work in these columns.
        # because array<str> deal with NA as "NA", the function can not recognize NA.
        obj = hl.if_else(mt.info[origFieldName][0] == 'NA',
                         hl.missing('float64'),
                         hl.float64(mt.info[origFieldName][0]))
        return obj
    
    # Called by annovar_vcf_preprocess
    def import_myvcf(self, path1):
        mt = hl.import_vcf(path1,reference_genome = 'GRCh37')
        return mt
    
    # Called by annovar_vcf_preprocess
    # filtering low complexity region
    def lcr_filter(self, mt, lcr_bed) :
        lcr = hl.import_bed(lcr_bed, reference_genome='GRCh37')
        mt = hl.filter_intervals(mt, lcr['interval'].collect(), keep=False)
        return mt
    
    # Called by annovar_vcf_preprocess
    # unkey 
    def mt_unkey(self, mt) :
        mt = mt.key_cols_by(); mt = mt.key_rows_by()
        return mt
    
    ### Filtering ###
    # Called by annovar_vcf_preprocess
    # Even though it is evident that a variant is het ,according to GT column,
    # some variant has not normalized PL
    # in other words, the variant does not have zero in het PL score.
    def filtering_myvcf(self, matrixtable, gnomadCutoff, GQCutoff, DPCutoff):
        mt = matrixtable    
        mt = mt.annotate_rows(gnomAD_r211_merged_all = hl.if_else(hl.is_defined(mt.info["gnomAD_r2.1.1_merged_all"]),
                                                          mt.info["gnomAD_r2.1.1_merged_all"], 0))
        mt = mt.filter_rows((mt.gnomAD_r211_merged_all <= gnomadCutoff) & # @arg : gnomadCutoff
                            (mt.info.MQ >= 40))
        mt = mt.annotate_entries(PLDIFF = "hom") # annotate PLDIFF
        mt = mt.annotate_entries(PLDIFFDP = "hom") # annotate PLDIFFDP

        mt = mt.filter_entries((mt.GQ > GQCutoff) & # @arg : GQCutoff
                               (mt.DP >= DPCutoff) & # @arg : DPCutoff
                               (mt.GT.is_hom_var()) & # recessive model
                               #(mt.PLDIFFDP >= PLDIFFDPCutoff) & # @arg : PLDIFFDPCutoff
                               (mt.PL[2] == 0), keep = True) # non normalized PL
        return mt
    
    def select_rows_from_mt(self, matrixtable) :
        mt = matrixtable
        mt = mt.select_rows(mt.locus,
                            mt.Ref,
                            mt.Alt,
                            mt.AC,
                            mt.AF,
                            mt.AN,
                            mt.gnomad_constraint_zscore,
                            mt.info["gnomAD_r2.1.1_merged_all"],
                            mt.info["gnomAD_r2.1_merged_hom"],
                            mt.info["MQ"],
                            mt.gnomad_genome_AF_eas,
                            mt.bravo_freeze8,
                            mt.oe_lof_gnomad_refGeneWithVer,
                            mt.oe_lof_lower_gnomad_refGeneWithVer,
                            mt.oe_lof_upper_gnomad_refGeneWithVer,
                            mt.oe_mis_lower_gnomad_refGeneWithVer,
                            mt.oe_mis_upper_gnomad_refGeneWithVer,
                            mt.pLi_exac_refGeneWithVer,
                            mt.pLi_gnomad_refGeneWithVer,
                            mt.MetaSVM_pred,
                            mt.MetaSVM_rankscore,
                            mt.MetaSVM_score,
                            mt.REVEL_rankscore,
                            mt.REVEL_score,
                            mt.spliceAI_DS_AG,
                            mt.spliceAI_DS_AL,
                            mt.spliceAI_DS_DG,
                            mt.spliceAI_DS_DL,
                            mt.AAChange,
                            mt.mis_z_gnomad_refGeneWithVer,
                            mt.mis_z_exac_refGeneWithVer,
                            mt.Func_refGeneWithVer,
                            mt.Gene_refGeneWithVer,
                            mt.ExonicFunc_refGeneWithVer,
                            mt.Gene_full_name_refGeneWithVer,
                            mt.CADD13,
                            mt.stopremoving_effect,
                            mt.stopremoving_kozak_strength,
                            mt.uAUGcreating_effect,
                            mt.uAUGcreating_kozak_strength,
                            mt.encode_CTCF_only,
                            mt.encode_DNase_H3K4me3,
                            mt.encode_PLS,
                            mt.encode_cCREs,
                            mt.encode_dELS,
                            mt.encode_pELS,
                            mt.ChromHMM_GM12878,
                            mt.ChromHMM_H1_hESC,
                            mt.ChromHMM_NHLF,
                            mt.vista_enhancers,
                            mt.fantom5_enhancers,
                            mt.jarvis)
        return mt
    
    # Called by annovar_vcf_preprocess
    def export_hailtable(self, matrixtable, path2): 
        mt = matrixtable
        ht = mt.entries(); ht = ht.key_by()
        ht.export(path2, header = True, delimiter = "\t")
    
    # main function 
    def annovar_vcf_preprocess(self, path1, lcr_bed, gnomadCutoff, GQCutoff, DPCutoff): 
        mt = self.import_myvcf(path1)
        mt = self.lcr_filter(mt, lcr_bed)
        mt = self.mt_unkey(mt)
        path2 = path1.replace(".vcf", ".rec_variant.table")
        gnomadCutoff = float(gnomadCutoff)
        GQCutoff = int(GQCutoff)
        DPCutoff = int(DPCutoff)
        mt = mt.annotate_rows(Ref = mt.alleles[0],
                              Alt = mt.alleles[1],
                              AC = self.float_float_preprocess(mt, "AC"),
                              AF = self.float_float_preprocess(mt, "AF"),
                              AN = hl.int32(mt.info.AN),
                              gnomad_constraint_zscore = self.str_float_preprocess(mt, "gnomad_constraint_zscore"),
                              gnomad_genome_AF_eas = self.str_float_preprocess(mt, "gnomad_genome_AF_eas"),
                              bravo_freeze8 = self.str_float_preprocess(mt, "bravo_freeze8"),
                              oe_lof_gnomad_refGeneWithVer = self.str_float_preprocess(mt, "oe_lof_gnomad.refGeneWithVer"),
                              oe_lof_lower_gnomad_refGeneWithVer = self.str_float_preprocess(mt, "oe_lof_lower_gnomad.refGeneWithVer"),
                              oe_lof_upper_gnomad_refGeneWithVer = self.str_float_preprocess(mt, "oe_lof_upper_gnomad.refGeneWithVer"),
                              oe_mis_lower_gnomad_refGeneWithVer = self.str_float_preprocess(mt, "oe_mis_lower_gnomad.refGeneWithVer"),
                              oe_mis_upper_gnomad_refGeneWithVer = self.str_float_preprocess(mt, "oe_mis_upper_gnomad.refGeneWithVer"),
                              pLi_exac_refGeneWithVer = self.str_float_preprocess(mt, "pLi_exac.refGeneWithVer"),
                              pLi_gnomad_refGeneWithVer = self.str_float_preprocess(mt, "pLi_gnomad.refGeneWithVer"),
                              MetaSVM_pred = self.str_str_preprocess(mt, "MetaSVM_pred"),
                              MetaSVM_rankscore = self.str_float_preprocess(mt, "MetaSVM_rankscore"),
                              MetaSVM_score = self.str_float_preprocess(mt, "MetaSVM_score"),
                              REVEL_rankscore = self.str_float_preprocess(mt, "REVEL_rankscore"),
                              REVEL_score = self.str_float_preprocess(mt, "REVEL_score"),
                              CADD13 = self.str_float_preprocess(mt, "CADD13_PHRED"),
                              spliceAI_DS_AG = self.str_float_preprocess(mt, "spliceAI_DS_AG"),
                              spliceAI_DS_AL = self.str_float_preprocess(mt, "spliceAI_DS_AL"),
                              spliceAI_DS_DG = self.str_float_preprocess(mt, "spliceAI_DS_DG"),
                              spliceAI_DS_DL = self.str_float_preprocess(mt, "spliceAI_DS_DL"),
                              csq = self.str_str_preprocess(mt, "csq"),
                              AAChange = self.str_str_preprocess(mt, "AAChange.refGeneWithVer"),
                              mis_z_gnomad_refGeneWithVer = self.str_float_preprocess(mt, "mis_z_gnomad.refGeneWithVer"),
                              mis_z_exac_refGeneWithVer = self.str_float_preprocess(mt, "mis_z_exac.refGeneWithVer"),
                              Func_refGeneWithVer = self.str_str_preprocess(mt, "Func.refGeneWithVer"),
                              Gene_refGeneWithVer = self.str_str_preprocess(mt, "Gene.refGeneWithVer"),
                              ExonicFunc_refGeneWithVer = self.str_str_preprocess(mt, "ExonicFunc.refGeneWithVer"),
                              Gene_full_name_refGeneWithVer = self.str_str_preprocess(mt, "Gene_full_name.refGeneWithVer"),
                              stopremoving_effect = self.str_str_preprocess(mt, "stopremoving_effect"),
                              stopremoving_kozak_strength = self.str_str_preprocess(mt, "stopremoving_kozak_strength"),
                              uAUGcreating_effect = self.str_str_preprocess(mt, "uAUGcreating_effect"),
                              uAUGcreating_kozak_strength = self.str_str_preprocess(mt, "uAUGcreating_kozak_strength"),
                              encode_CTCF_only = self.str_str_preprocess(mt, "encode_CTCF_only"),
                              encode_DNase_H3K4me3 = self.str_str_preprocess(mt, "encode_DNase_H3K4me3"),
                              encode_PLS = self.str_str_preprocess(mt, "encode_PLS"),
                              encode_cCREs = self.str_str_preprocess(mt, "encode_cCREs"),
                              encode_dELS = self.str_str_preprocess(mt, "encode_dELS"),
                              encode_pELS = self.str_str_preprocess(mt, "encode_pELS"),
                              ChromHMM_GM12878 = self.str_str_preprocess(mt, "ChromHMM_GM12878"),
                              ChromHMM_H1_hESC = self.str_str_preprocess(mt, "ChromHMM_H1_hESC"),
                              ChromHMM_NHLF = self.str_str_preprocess(mt, "ChromHMM_NHLF"),
                              vista_enhancers = self.str_str_preprocess(mt, "vista_enhancers"),
                              fantom5_enhancers = self.str_str_preprocess(mt, "fantom5_enhancers"),
                              jarvis = self.str_float_preprocess(mt, "jarvis"))
        
        mt = self.filtering_myvcf(mt, gnomadCutoff, GQCutoff, DPCutoff)
        mt = self.select_rows_from_mt(mt)
        self.export_hailtable(mt, path2)

myobj = AnnovarWithHail()
myobj.annovar_vcf_preprocess(input_file, lcr_bed, MAF, 10, 4)
hl.stop()
