#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 00:20:57 2023

@author: ptruong
"""

import re
s = "{0.01929014233974128}" # the original string
r = r"\{(\d*\.?\d*)\}" # the regular expression to match a number inside curly brackets
t = "[\\1]" # the replacement string to replace the matched number with square brackets
result = re.sub(r, t, s) # the result of the replacement
print(result) # print the result

str = """{
    "model":{
      "EFV":    {
   {0.01929014233974128} 
       {0.01622765628659735} 
       {0.01769443322434847} 
       {0.01623235001696338} 
       {0.01622972861556908} 
       {0.01365311115696728} 
       {0.01488718145152558} 
       {0.0136570602190682} 
       {0.01769333756976861} 
       {0.01488435883925134} 
       {0.01622971850752687} 
       {0.01488866403069913} 
       {0.01623137334257151} 
       {0.01365449476855555} 
       {0.01488869012427688} 
       {0.01365844423085596} 
       {0.01929347190999687}
       },
       "Equilibrium frequency estimator":"Corrected 3x4 frequency estimator",
       "ID":"simulator.substitution_model",
       "Q":    {
   {"", "simulator.substitution_model.theta_AC*beta*0.23367778810339", "simulator.substitution_model.theta_AG*alpha*0.2547993342097017", "simulator.substitution_model.theta_AT*beta*0.2337453776869083", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""} 
       {"simulator.substitution_model.theta_AC*beta*0.2777775", "", "simulator.substitution_model.theta_CG*beta*0.2547993342097017", "simulator.substitution_model.theta_CT*alpha*0.2337453776869083", "", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2990038356551464", "", "", "", "", "", "", "", "", "", "", "", ""}
       },
       "alphabet":    {
   {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"} 
       },
      "bases":    {
   {"A", "C", "G", "T"} 
       },
       "branch-length-string":"0.05888392960524094*simulator.substitution_model.theta_GT*alpha+0.2916374429049448*simulator.substitution_model.theta_GT*beta+0.1354994125280253*simulator.substitution_model.theta_CT*alpha+0.2151052928084076*simulator.substitution_model.theta_CT*beta+0.05886675304670659*simulator.substitution_model.theta_CG*alpha+0.2916505948090215*simulator.substitution_model.theta_CG*beta+0.0717822104730694*simulator.substitution_model.theta_AT*alpha+0.2786537842410506*simulator.substitution_model.theta_AT*beta+0.1190703566772678*simulator.substitution_model.theta_AG*alpha+0.2524309597208058*simulator.substitution_model.theta_AG*beta+0.0876166671369822*simulator.substitution_model.theta_AC*alpha+0.262815720090049*simulator.substitution_model.theta_AC*beta",
      "canonical":0,
      "constrain-branch-length":"models.codon.MG_REV.ConstrainBranchLength",
      "defineQ":"models.codon.MG_REV._DefineQ",
      "description":"The Muse-Gaut 94 codon-substitution model coupled with the general time reversible (GTR) model of nucleotide substitution",
      "efv-id":"simulator.substitution_model_pi",
      "frequency-estimator":"frequencies.empirical.corrected.CF3x4",
      "get-branch-length":"",
      "matrix-id":"simulator.substitution_model_Q",
      "parameters":{
        "empirical":9,
        "global":{
          "Substitution rate from nucleotide A to nucleotide C":"simulator.substitution_model.theta_AC",
          "Substitution rate from nucleotide A to nucleotide G":"simulator.substitution_model.theta_AG",
          "Substitution rate from nucleotide A to nucleotide T":"simulator.substitution_model.theta_AT",
          "Substitution rate from nucleotide C to nucleotide G":"simulator.substitution_model.theta_CG",
          "Substitution rate from nucleotide C to nucleotide T":"simulator.substitution_model.theta_CT",
          "Substitution rate from nucleotide G to nucleotide T":"simulator.substitution_model.theta_GT"
         },
        "local":{
          "non-synonymous rate":"beta",
          "synonymous rate":"alpha"
         }
       },
       "post-definition":"models.generic.post.definition",
      "q_ij":"models.codon.MG_REV._GenerateRate",
      "reversible":1,
      "set-branch-length":"models.codon.MG_REV.set_branch_length",
      "stop":    {
   {"TAA", "TAG", "TGA"} 
       },
       "time":"models.DNA.generic.Time",
      "translation-table":{
        "AAA":"K",
        "AAC":"N",
        "AAG":"K",
        "AAT":"N",
        "ACA":"T",
        "ACC":"T",
        "ACG":"T",
        "ACT":"T",
        "AGA":"R",
        "AGC":"S",
        "AGG":"R",
        "AGT":"S",
        "ATA":"I",
        "ATC":"I",
        "ATG":"M",
        "ATT":"I",
        "CAA":"Q",
        "CAC":"H",
        "CAG":"Q",
        "CAT":"H",
        "CCA":"P",
        "CCC":"P",
        "CCG":"P",
        "CCT":"P",
        "CGA":"R",
        "CGC":"R",
        "CGG":"R",
        "CGT":"R",
        "CTA":"L",
        "CTC":"L",
        "CTG":"L",
        "CTT":"L",
        "GAA":"E",
        "GAC":"D",
        "GAG":"E",
        "GAT":"D",
        "GCA":"A",
        "GCC":"A",
        "GCG":"A",
        "GCT":"A",
        "GGA":"G",
        "GGC":"G",
        "GGG":"G",
        "GGT":"G",
        "GTA":"V",
        "GTC":"V",
        "GTG":"V",
        "GTT":"V",
        "TAA":"X",
        "TAC":"Y",
        "TAG":"X",
        "TAT":"Y",
        "TCA":"S",
        "TCC":"S",
        "TCG":"S",
        "TCT":"S",
        "TGA":"X",
        "TGC":"C",
        "TGG":"W",
        "TGT":"C",
        "TTA":"L",
        "TTC":"F",
        "TTG":"L",
        "TTT":"F"
       },
       "type":"local"
     },
     "simulator.site.profile":{
      "rates":    {
   {0.77, 0, 0, 0, 0} 
       {2.98, 0, 0.4, 0, 0.14} 
       {1.46, 1.36, 1.36, 1.36, 1.36} 
       {0.51, 0, 0, 0, 0} 
       {1.34, 0.03, 0.52, 0, 0.02} 
       {0.05, 0.17, 0.17, 0.17, 0.17} 
       {4.24, 0, 0, 0, 0} 
       {0.5600000000000001, 1.12, 1.12, 1.12, 1.12} 
       {0.92, 0, 0, 0, 0} 
       {0.6899999999999999, 0.5600000000000001, 0, 0.43, 0.08} 
       {0.42, 0.09, 0.09, 0.09, 0.09} 
       {1.81, 2.2, 2.2, 2.2, 2.2} 
       {0.32, 0.6, 0.6, 0.6, 0.6} 
       {1.33, 0, 0, 0, 0} 
       {0.2, 0, 0, 0, 0} 
       {0.3, 0, 0, 0, 0} 
       {0.68, 0.02, 0.02, 0.02, 0.02} 
       {0.86, 0.04, 0, 1.84, 0.05} 
       {0.63, 0, 0, 0.29, 0} 
       },
       "variables":{
        "0":"simulator.omega.class0",
        "1":"simulator.omega.class1",
        "2":"simulator.omega.class2",
        "3":"simulator.omega.class3"
       }
     },
     "tree":{
      "annotated_string":"((((((N4{GROUP_3},N5{GROUP_3})N69{GROUP_2},(N6{GROUP_2},N7{GROUP_0})N70{GROUP_3})N101{GROUP_3},((N8{GROUP_2},N9{GROUP_3})N71{GROUP_3},(N10{GROUP_3},N11{GROUP_3})N72{GROUP_2})N102{GROUP_2})N117{GROUP_2},(((N12{GROUP_1},N13{GROUP_1})N73{GROUP_1},(N14{GROUP_1},N15{GROUP_1})N74{GROUP_1})N103{GROUP_0},((N16{GROUP_0},N17{GROUP_0})N75{GROUP_0},(N18{GROUP_0},N19{GROUP_1})N76{GROUP_0})N104{GROUP_0})N118{GROUP_1})N125{GROUP_2},((((N20{GROUP_1},N21{GROUP_1})N77{GROUP_1},(N22{GROUP_0},N23{GROUP_0})N78{GROUP_0})N105{GROUP_0},((N24{GROUP_1},N25{GROUP_0})N79{GROUP_0},(N26{GROUP_0},N27{GROUP_0})N80{GROUP_3})N106{GROUP_3})N119{GROUP_2},(((N28{GROUP_1},N29{GROUP_3})N81{GROUP_3},(N30{GROUP_3},N31{GROUP_0})N82{GROUP_3})N107{GROUP_3},((N32{GROUP_0},N33{GROUP_0})N83{GROUP_2},(N34{GROUP_3},N35{GROUP_2})N84{GROUP_1})N108{GROUP_3})N120{GROUP_3})N126{GROUP_0})N129{GROUP_0},(((((N36{GROUP_0},N37{GROUP_1})N85{GROUP_1},(N38{GROUP_2},N39{GROUP_1})N86{GROUP_2})N109{GROUP_2},((N40{GROUP_0},N41{GROUP_2})N87{GROUP_1},(N42{GROUP_1},N43{GROUP_1})N88{GROUP_2})N110{GROUP_0})N121{GROUP_2},(((N44{GROUP_2},N45{GROUP_0})N89{GROUP_0},(N46{GROUP_0},N47{GROUP_3})N90{GROUP_1})N111{GROUP_2},((N48{GROUP_2},N49{GROUP_3})N91{GROUP_1},(N50{GROUP_0},N51{GROUP_0})N92{GROUP_0})N112{GROUP_0})N122{GROUP_0})N127{GROUP_0},((((N52{GROUP_3},N53{GROUP_3})N93{GROUP_1},(N54{GROUP_3},N55{GROUP_3})N94{GROUP_1})N113{GROUP_0},((N56{GROUP_0},N57{GROUP_2})N95{GROUP_2},(N58{GROUP_2},N59{GROUP_0})N96{GROUP_0})N114{GROUP_3})N123{GROUP_3},(((N60{GROUP_1},N61{GROUP_1})N97{GROUP_0},(N62{GROUP_1},N63{GROUP_3})N98{GROUP_1})N115{GROUP_2},((N64{GROUP_2},N65{GROUP_3})N99{GROUP_2},((N0{GROUP_1},N1{GROUP_2})N67{GROUP_3},(N2{GROUP_1},N3{GROUP_1})N68{GROUP_1})N100{GROUP_3})N116{GROUP_2})N124{GROUP_1})N128{GROUP_1})N130{GROUP_0},N66{GROUP_0});",
      "branch length":{
        "N0":0.09683796328447575,
        "N1":0.176288090497566,
        "N10":0.08924779452315575,
        "N100":0.009420853927403763,
        "N101":0.2651641240293992,
        "N102":0.148155084665296,
        "N103":0.02339025087803294,
        "N104":0.09537525969712085,
        "N105":0.1597008906424866,
        "N106":0.1282248698858558,
        "N107":0.04509091951102472,
        "N108":0.07812869165398892,
        "N109":0.003839970586509975,
        "N11":0.02554457109259972,
        "N110":0.06744838995912472,
        "N111":0.004671638351592623,
        "N112":0.1472213934081857,
        "N113":0.02817948195486256,
        "N114":0.1145883863445523,
        "N115":0.03569201124102021,
        "N116":0.1179067920348875,
        "N117":0.05189725576830864,
        "N118":0.08248101929435084,
        "N119":0.04018184535324473,
        "N12":0.008610665160227935,
        "N120":0.06443056944800284,
        "N121":0.01343706731661185,
        "N122":0.165041829485691,
        "N123":0.04709209981921693,
        "N124":0.2461586447709646,
        "N125":0.04517093684168755,
        "N126":0.01806521728455023,
        "N127":0.3678363624003491,
        "N128":0.2801592257422743,
        "N129":0.1473676069887016,
        "N13":0.06148056069761299,
        "N130":0.02246580224145033,
        "N14":0.1307313489917125,
        "N15":0.0424873492734112,
        "N16":0.2716823769646388,
        "N17":0.04824487440516276,
        "N18":0.1073793780269826,
        "N19":0.02512749368725543,
        "N2":0.04550739881556531,
        "N20":0.07901647780237897,
        "N21":0.08467138130509273,
        "N22":0.002055209168422354,
        "N23":0.1030125188924388,
        "N24":0.07139749640544509,
        "N25":0.01083730482432979,
        "N26":0.109940318093655,
        "N27":0.01537831276625502,
        "N28":0.1042028295336127,
        "N29":0.06879534879583676,
        "N3":0.02201800229780494,
        "N30":0.0281719652003367,
        "N31":0.1714194701332908,
        "N32":0.3981312708731555,
        "N33":0.02161944237806132,
        "N34":0.08060266107451972,
        "N35":0.02901566817242126,
        "N36":0.1829710480484212,
        "N37":0.01496601538509164,
        "N38":0.003822133510925323,
        "N39":0.1833785409047971,
        "N4":0.320482587418795,
        "N40":0.2279343107384088,
        "N41":0.01169833106368832,
        "N42":0.1502680379300581,
        "N43":0.1798956981591125,
        "N44":0.3146892906370151,
        "N45":0.1665904772970213,
        "N46":0.27011517108781,
        "N47":0.02541530075905582,
        "N48":0.08580783225492532,
        "N49":0.3407116738896538,
        "N5":0.0767934986309215,
        "N50":0.0753353232522052,
        "N51":0.02994516664648106,
        "N52":0.05797517093606119,
        "N53":0.108697848747349,
        "N54":0.06978552552636991,
        "N55":0.007491099548399317,
        "N56":0.03098264634967597,
        "N57":0.004437471636468425,
        "N58":0.01715543967488304,
        "N59":0.02446898973774961,
        "N6":0.0579015592298889,
        "N60":0.07680302238726837,
        "N61":0.05235498807694206,
        "N62":0.2796277868205672,
        "N63":0.03136952720325668,
        "N64":0.1166301438710106,
        "N65":0.1246901618223763,
        "N66":0.08695158227067276,
        "N67":0.03347657442701567,
        "N68":0.0301044092138499,
        "N69":0.004734476716270237,
        "N7":0.08865935157270551,
        "N70":0.2560093020590626,
        "N71":0.08304650005021595,
        "N72":0.04680148998689951,
        "N73":0.2198322097819825,
        "N74":0.1181823147012625,
        "N75":0.007630052239767903,
        "N76":0.008346126846801654,
        "N77":0.065738710164002,
        "N78":0.000839607030336947,
        "N79":0.1913812065541395,
        "N8":0.1027485044850314,
        "N80":0.1605316442026653,
        "N81":0.04524344596268309,
        "N82":0.01859878266984774,
        "N83":0.3887066295861277,
        "N84":0.1082902226880476,
        "N85":0.04577061720267689,
        "N86":0.2830561375294785,
        "N87":0.04588579711374348,
        "N88":0.07602961662461,
        "N89":0.03003543233269226,
        "N9":0.05577636323342411,
        "N90":0.1034007514058179,
        "N91":0.01154206186888107,
        "N92":0.04444395504489117,
        "N93":0.09969689289010175,
        "N94":0.02619872991300569,
        "N95":0.07541900464536098,
        "N96":0.2962896870110845,
        "N97":0.04962236094266016,
        "N98":0.08004586215727705,
        "N99":0.0009703235398180994
       },
      "file":"/Users/sergei/Dropbox/Work/FEL-contrast/datasets/omnibus-multi/sims.0.nwk",
      "meta":{
       },
      "model_list":    {
   {"GROUP_0", "GROUP_1", "GROUP_2", "GROUP_3"} 
       },
      "model_map":{
        "N0":"GROUP_1",
        "N1":"GROUP_2",
        "N10":"GROUP_3",
        "N100":"GROUP_3",
        "N101":"GROUP_3",
        "N102":"GROUP_2",
        "N103":"GROUP_0",
        "N104":"GROUP_0",
        "N105":"GROUP_0",
        "N106":"GROUP_3",
        "N107":"GROUP_3",
        "N108":"GROUP_3",
        "N109":"GROUP_2",
        "N11":"GROUP_3",
        "N110":"GROUP_0",
        "N111":"GROUP_2",
        "N112":"GROUP_0",
        "N113":"GROUP_0",
        "N114":"GROUP_3",
        "N115":"GROUP_2",
        "N116":"GROUP_2",
        "N117":"GROUP_2",
        "N118":"GROUP_1",
        "N119":"GROUP_2",
       },
      "partitioned":{
        "N0":"leaf",
        "N1":"leaf",
        "N10":"leaf",
        "N100":"internal",
        "N101":"internal",
        "N102":"internal",
        "N103":"internal",
        "N104":"internal",
        "N105":"internal",
        "N106":"internal",
        "N107":"internal",
        "N108":"internal",
        "N109":"internal",
        "N11":"leaf",
        "N110":"internal",
        "N111":"internal",
        "N112":"internal",
        "N113":"internal",
        "N114":"internal",
        "N115":"internal",
        "N116":"internal",
        "N117":"internal",
        "N118":"internal",
        "N119":"internal",
       },
      "rooted":0,
      "string":"((((((N4,N5)N69,(N6,N7)N70)N101,((N8,N9)N71,(N10,N11)N72)N102)N117,(((N12,N13)N73,(N14,N15)N74)N103,((N16,N17)N75,(N18,N19)N76)N104)N118)N125,((((N20,N21)N77,(N22,N23)N78)N105,((N24,N25)N79,(N26,N27)N80)N106)N119,(((N28,N29)N81,(N30,N31)N82)N107,((N32,N33)N83,(N34,N35)N84)N108)N120)N126)N129,(((((N36,N37)N85,(N38,N39)N86)N109,((N40,N41)N87,(N42,N43)N88)N110)N121,(((N44,N45)N89,(N46,N47)N90)N111,((N48,N49)N91,(N50,N51)N92)N112)N122)N127,((((N52,N53)N93,(N54,N55)N94)N113,((N56,N57)N95,(N58,N59)N96)N114)N123,(((N60,N61)N97,(N62,N63)N98)N115,((N64,N65)N99,((N0,N1)N67,(N2,N3)N68)N100)N116)N124)N128)N130,N66)",   "string_with_lengths":"((((((N4:0.320482587418795,N5:0.0767934986309215)N69:0.004734476716270237,(N6:0.0579015592298889,N7:0.08865935157270551)N70:0.2560093020590626)N101:0.2651641240293992,((N8:0.1027485044850314,N9:0.05577636323342411)N71:0.08304650005021595,(N10:0.08924779452315575,N11:0.02554457109259972)N72:0.04680148998689951)N102:0.148155084665296)N117:0.05189725576830864,(((N12:0.008610665160227935,N13:0.06148056069761299)N73:0.2198322097819825,(N14:0.1307313489917125,N15:0.0424873492734112)N74:0.1181823147012625)N103:0.02339025087803294,((N16:0.2716823769646388,N17:0.04824487440516276)N75:0.007630052239767903,(N18:0.1073793780269826,N19:0.02512749368725543)N76:0.008346126846801654)N104:0.09537525969712085)N118:0.08248101929435084)N125:0.04517093684168755,((((N20:0.07901647780237897,N21:0.08467138130509273)N77:0.065738710164002,(N22:0.002055209168422354,N23:0.1030125188924388)N78:0.000839607030336947)N105:0.1597008906424866,((N24:0.07139749640544509,N25:0.01083730482432979)N79:0.1913812065541395,(N26:0.109940318093655,N27:0.01537831276625502)N80:0.1605316442026653)N106:0.1282248698858558)N119:0.04018184535324473,(((N28:0.1042028295336127,N29:0.06879534879583676)N81:0.04524344596268309,(N30:0.0281719652003367,N31:0.1714194701332908)N82:0.01859878266984774)N107:0.04509091951102472,((N32:0.3981312708731555,N33:0.02161944237806132)N83:0.3887066295861277,(N34:0.08060266107451972,N35:0.02901566817242126)N84:0.1082902226880476)N108:0.07812869165398892)N120:0.06443056944800284)N126:0.01806521728455023)N129:0.1473676069887016,(((((N36:0.1829710480484212,N37:0.01496601538509164)N85:0.04577061720267689,(N38:0.003822133510925323,N39:0.1833785409047971)N86:0.2830561375294785)N109:0.003839970586509975,((N40:0.2279343107384088,N41:0.01169833106368832)N87:0.04588579711374348,(N42:0.1502680379300581,N43:0.1798956981591125)N88:0.07602961662461)N110:0.06744838995912472)N121:0.01343706731661185,(((N44:0.3146892906370151,N45:0.1665904772970213)N89:0.03003543233269226,(N46:0.27011517108781,N47:0.02541530075905582)N90:0.1034007514058179)N111:0.004671638351592623,((N48:0.08580783225492532,N49:0.3407116738896538)N91:0.01154206186888107,(N50:0.0753353232522052,N51:0.02994516664648106)N92:0.04444395504489117)N112:0.1472213934081857)N122:0.165041829485691)N127:0.3678363624003491,((((N52:0.05797517093606119,N53:0.108697848747349)N93:0.09969689289010175,(N54:0.06978552552636991,N55:0.007491099548399317)N94:0.02619872991300569)N113:0.02817948195486256,((N56:0.03098264634967597,N57:0.004437471636468425)N95:0.07541900464536098,(N58:0.01715543967488304,N59:0.02446898973774961)N96:0.2962896870110845)N114:0.1145883863445523)N123:0.04709209981921693,(((N60:0.07680302238726837,N61:0.05235498807694206)N97:0.04962236094266016,(N62:0.2796277868205672,N63:0.03136952720325668)N98:0.08004586215727705)N115:0.03569201124102021,((N64:0.1166301438710106,N65:0.1246901618223763)N99:0.0009703235398180994,((N0:0.09683796328447575,N1:0.176288090497566)N67:0.03347657442701567,(N2:0.04550739881556531,N3:0.02201800229780494)N68:0.0301044092138499)N100:0.009420853927403763)N116:0.1179067920348875)N124:0.2461586447709646)N128:0.2801592257422743)N130:0.02246580224145033,N66:0.08695158227067276)"
     }
   }
       """


# Open the file for reading
with open('settings_sample.file', 'r') as file:
    # Read the entire file contents into a string
    file_contents = file.read()
    
str = file_contents 

def replace_single_value_curly_brackets(s):
    r = r"\{(\d*\.?\d*)\}" # the regular expression to match a number inside curly brackets
    t = "[\\1]" # the replacement string to replace the matched number with square brackets
    result = re.sub(r, t, s) # the result of the replacement
    return result # print the result
     
def replace_vector_value_curly_brackets(s):
    #s = "{0.77, 0, 0, 0, 0}" # the original string  
    r = r"\{([\d\., ]*)\}" # the regular expression to match a comma-separated list of numbers inside curly brackets
    t = "[\\1]" # the replacement string to replace the matched list with square brackets
    result = re.sub(r, t, s) # the result of the replacement
    return result # print the result

def fix_Q_key_curly_brackets(s):
    # Define a regular expression to match and replace inner curly brackets
    pattern = r'(\{.*?\})'
    # Replace inner curly brackets with square brackets
    output_text = re.sub(pattern, lambda m: m.group(0).replace('{', '[').replace('}', ']'), s)
    return output_text




def fix_EFV_list(s):
    r = r"\n\s+" # the regular expression to match a newline character followed by one or more whitespace characters
    t = ", " # the replacement string to replace the matched part with a comma and a space
    result = re.sub(r, t, s) # the result of the replacement
    return result # print the result

res = replace_single_value_curly_brackets(str)
res = replace_vector_value_curly_brackets(res)
#res = fix_EFV_list(res)
res = fix_Q_key_curly_brackets(res)

# this bugs when we are working with fix_Q_key_curly_brackets()
res = re.sub(r'(?<=\])\s+(?=\[)', ',', res)




# Define a regular expression pattern to match elements within square brackets
#pattern = r'\[([^]]*)\]'

# Replace elements within square brackets with elements followed by a comma
#output_text = re.sub(pattern, r'[\1],', res)
#print(res)


# Use regex to replace { and } with [ and ]
#res = re.sub(r'\{', '[', str)
#res = re.sub(r'\}', ']', res)

# Print the modified string
print(res)

smaller_sample = """
{
    "model":{
      "EFV":    {
   [0.01929014233974128],[0.01622765628659735],[0.01769443322434847],[0.01623235001696338],[0.01622972861556908],[0.01365311115696728],[0.01488718145152558],[0.0136570602190682],[0.01769333756976861],[0.01488435883925134],[0.01622971850752687],[0.01488866403069913],[0.01623137334257151],[0.01365449476855555],[0.01488869012427688],[0.01365844423085596],[0.01929347190999687]
       },
       "Equilibrium frequency estimator":"Corrected 3x4 frequency estimator",
       "ID":"simulator.substitution_model",
       "Q":    {
   ["", "simulator.substitution_model.theta_AC*beta*0.23367778810339", "simulator.substitution_model.theta_AG*alpha*0.2547993342097017", "simulator.substitution_model.theta_AT*beta*0.2337453776869083", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""],["simulator.substitution_model.theta_AC*beta*0.2777775", "", "simulator.substitution_model.theta_CG*beta*0.2547993342097017", "simulator.substitution_model.theta_CT*alpha*0.2337453776869083", "", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2990038356551464", "", "", "", "", "", "", "", "", "", "", "", ""]
       }}}
"""

# Define a list of keys to transform
keys_to_transform = ["EFV", "Q", "alphabet", "bases","stop", 
                     "rates", "meta", "model_list"]

# Use a loop to transform the keys in the sample
transformed_sample = res

for key in keys_to_transform:
    # Use regex to perform the transformation for the specific key
    pattern = rf'"{key}"\s*:\s*{{([^}}]*)}}'
    replacement = rf'"{key}": [\1]'
    transformed_sample = re.sub(pattern, replacement, transformed_sample)


transformed_sample = re.sub(r',(?=\s*})', '', transformed_sample)
print(transformed_sample)

import json

json.loads(transformed_sample)