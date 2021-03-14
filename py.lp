\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
 [ 60 state[0] ^2 + 60 state[1] ^2 + 60 state[2] ^2 + 60 state[3] ^2
   + 60 state[4] ^2 + 60 state[5] ^2 + 60 state[6] ^2 + 60 state[7] ^2
   + 60 state[8] ^2 + 60 state[9] ^2 + 60 state[10] ^2 + 60 state[11] ^2
   + 60 state[12] ^2 + 60 state[13] ^2 + 60 state[14] ^2 + 60 state[15] ^2
   + 60 state[16] ^2 + 60 state[17] ^2 + 60 state[18] ^2
   + 12000 state[19] ^2 + 60 state[20] ^2 + 60 state[21] ^2
   + 60 state[22] ^2 + 60 state[23] ^2 + 60 state[24] ^2 + 60 state[25] ^2
   + 60 state[26] ^2 + 60 state[27] ^2 + 60 state[28] ^2 + 60 state[29] ^2
   + 60 state[30] ^2 + 60 state[31] ^2 + 60 state[32] ^2 + 60 state[33] ^2
   + 60 state[34] ^2 + 60 state[35] ^2 + 60 state[36] ^2 + 60 state[37] ^2
   + 60 state[38] ^2 + 12000 state[39] ^2 + 2 state[40] ^2 + 2 state[41] ^2
   + 2 state[42] ^2 + 2 state[43] ^2 + 2 state[44] ^2 + 2 state[45] ^2
   + 2 state[46] ^2 + 2 state[47] ^2 + 2 state[48] ^2 + 2 state[49] ^2
   + 2 state[50] ^2 + 2 state[51] ^2 + 2 state[52] ^2 + 2 state[53] ^2
   + 2 state[54] ^2 + 2 state[55] ^2 + 2 state[56] ^2 + 2 state[57] ^2
   + 2 state[58] ^2 + 400 state[59] ^2 + control[0] ^2 + control[1] ^2
   + control[2] ^2 + control[3] ^2 + control[4] ^2 + control[5] ^2
   + control[6] ^2 + control[7] ^2 + control[8] ^2 + control[9] ^2
   + control[10] ^2 + control[11] ^2 + control[12] ^2 + control[13] ^2
   + control[14] ^2 + control[15] ^2 + control[16] ^2 + control[17] ^2
   + control[18] ^2 + control[19] ^2 + control[20] ^2 + control[21] ^2
   + control[22] ^2 + control[23] ^2 + control[24] ^2 + control[25] ^2
   + control[26] ^2 + control[27] ^2 + control[28] ^2 + control[29] ^2
   + control[30] ^2 + control[31] ^2 + control[32] ^2 + control[33] ^2
   + control[34] ^2 + control[35] ^2 + control[36] ^2 + control[37] ^2
   + 0.01 control[38] ^2 + 0.01 control[39] ^2 + 0.01 control[40] ^2
   + 0.01 control[41] ^2 + 0.01 control[42] ^2 + 0.01 control[43] ^2
   + 0.01 control[44] ^2 + 0.01 control[45] ^2 + 0.01 control[46] ^2
   + 0.01 control[47] ^2 + 0.01 control[48] ^2 + 0.01 control[49] ^2
   + 0.01 control[50] ^2 + 0.01 control[51] ^2 + 0.01 control[52] ^2
   + 0.01 control[53] ^2 + 0.01 control[54] ^2 + 0.01 control[55] ^2
   + 0.01 control[56] ^2 ] / 2 
Subject To
 R0: - state[0] + state[1] - 0.00261806778230845 control[0] = 0
 R1: - state[20] + state[21] - 7.8542033469253367e-04 state[40]
   - 0.00261806778230845 control[19] = 0
 R2: - state[40] + state[41] + 0.0104166666666667 state[60]
   + 0.0347222222222222 control[19] = 0
 R3: - state[60] + state[61] - 0.05 control[38] = 0
 R4: 0.0523613556461689 control[0] >= -0.0157084066938507
 R5: 0.009375 state[60] + 0.0836113556461689 control[19]
   + 0.045 control[38] >= -0.3
 R6: 0.0523613556461689 control[0] <= 0.2842915933061493
 R7: 0.009375 state[60] + 0.0836113556461689 control[19]
   + 0.045 control[38] <= 0.3
 R8: - 0.6 control[0] + control[19] <= 0.18
 R9: 0.6 control[0] + control[19] >= -0.18
 R10: control[38] - mode[0] >= -1
 R11: control[38] + mode[0] <= 1
 R12: control[38] - 5 mode[19] >= -5
 R13: 0.6 control[0] - control[19] - mode[19] >= -1.18
 R14: 0.6 control[0] - control[19] + mode[19] <= 0.82
 R15: control[38] + 5 mode[38] <= 5
 R16: 0.6 control[0] + control[19] - mode[38] >= -1.18
 R17: 0.6 control[0] + control[19] + mode[38] <= 0.82
 R18: mode[0] + mode[19] + mode[38] = 1
 R19: - state[1] + state[2] - 0.00261806778230845 control[1] = 0
 R20: - state[21] + state[22] - 7.8542033469253367e-04 state[41]
   - 0.00261806778230845 control[20] = 0
 R21: - state[41] + state[42] + 0.0104166666666667 state[61]
   + 0.0347222222222222 control[20] = 0
 R22: - state[61] + state[62] - 0.05 control[39] = 0
 R23: 0.0523613556461689 control[1] >= -0.0157084066938507
 R24: 0.009375 state[61] + 0.0836113556461689 control[20]
   + 0.045 control[39] >= -0.3
 R25: 0.0523613556461689 control[1] <= 0.2842915933061493
 R26: 0.009375 state[61] + 0.0836113556461689 control[20]
   + 0.045 control[39] <= 0.3
 R27: - 0.6 control[1] + control[20] <= 0.18
 R28: 0.6 control[1] + control[20] >= -0.18
 R29: control[39] - mode[1] >= -1
 R30: control[39] + mode[1] <= 1
 R31: control[39] - 5 mode[20] >= -5
 R32: 0.6 control[1] - control[20] - mode[20] >= -1.18
 R33: 0.6 control[1] - control[20] + mode[20] <= 0.82
 R34: control[39] + 5 mode[39] <= 5
 R35: 0.6 control[1] + control[20] - mode[39] >= -1.18
 R36: 0.6 control[1] + control[20] + mode[39] <= 0.82
 R37: mode[1] + mode[20] + mode[39] = 1
 R38: - state[2] + state[3] - 0.00261806778230845 control[2] = 0
 R39: - state[22] + state[23] - 7.8542033469253367e-04 state[42]
   - 0.00261806778230845 control[21] = 0
 R40: - state[42] + state[43] + 0.0104166666666667 state[62]
   + 0.0347222222222222 control[21] = 0
 R41: - state[62] + state[63] - 0.05 control[40] = 0
 R42: 0.0523613556461689 control[2] >= -0.0157084066938507
 R43: 0.009375 state[62] + 0.0836113556461689 control[21]
   + 0.045 control[40] >= -0.3
 R44: 0.0523613556461689 control[2] <= 0.2842915933061493
 R45: 0.009375 state[62] + 0.0836113556461689 control[21]
   + 0.045 control[40] <= 0.3
 R46: - 0.6 control[2] + control[21] <= 0.18
 R47: 0.6 control[2] + control[21] >= -0.18
 R48: control[40] - mode[2] >= -1
 R49: control[40] + mode[2] <= 1
 R50: control[40] - 5 mode[21] >= -5
 R51: 0.6 control[2] - control[21] - mode[21] >= -1.18
 R52: 0.6 control[2] - control[21] + mode[21] <= 0.82
 R53: control[40] + 5 mode[40] <= 5
 R54: 0.6 control[2] + control[21] - mode[40] >= -1.18
 R55: 0.6 control[2] + control[21] + mode[40] <= 0.82
 R56: mode[2] + mode[21] + mode[40] = 1
 R57: - state[3] + state[4] - 0.00261806778230845 control[3] = 0
 R58: - state[23] + state[24] - 7.8542033469253367e-04 state[43]
   - 0.00261806778230845 control[22] = 0
 R59: - state[43] + state[44] + 0.0104166666666667 state[63]
   + 0.0347222222222222 control[22] = 0
 R60: - state[63] + state[64] - 0.05 control[41] = 0
 R61: 0.0523613556461689 control[3] >= -0.0157084066938507
 R62: 0.009375 state[63] + 0.0836113556461689 control[22]
   + 0.045 control[41] >= -0.3
 R63: 0.0523613556461689 control[3] <= 0.2842915933061493
 R64: 0.009375 state[63] + 0.0836113556461689 control[22]
   + 0.045 control[41] <= 0.3
 R65: - 0.6 control[3] + control[22] <= 0.18
 R66: 0.6 control[3] + control[22] >= -0.18
 R67: control[41] - mode[3] >= -1
 R68: control[41] + mode[3] <= 1
 R69: control[41] - 5 mode[22] >= -5
 R70: 0.6 control[3] - control[22] - mode[22] >= -1.18
 R71: 0.6 control[3] - control[22] + mode[22] <= 0.82
 R72: control[41] + 5 mode[41] <= 5
 R73: 0.6 control[3] + control[22] - mode[41] >= -1.18
 R74: 0.6 control[3] + control[22] + mode[41] <= 0.82
 R75: mode[3] + mode[22] + mode[41] = 1
 R76: - state[4] + state[5] - 0.00261806778230845 control[4] = 0
 R77: - state[24] + state[25] - 7.8542033469253367e-04 state[44]
   - 0.00261806778230845 control[23] = 0
 R78: - state[44] + state[45] + 0.0104166666666667 state[64]
   + 0.0347222222222222 control[23] = 0
 R79: - state[64] + state[65] - 0.05 control[42] = 0
 R80: 0.0523613556461689 control[4] >= -0.0157084066938507
 R81: 0.009375 state[64] + 0.0836113556461689 control[23]
   + 0.045 control[42] >= -0.3
 R82: 0.0523613556461689 control[4] <= 0.2842915933061493
 R83: 0.009375 state[64] + 0.0836113556461689 control[23]
   + 0.045 control[42] <= 0.3
 R84: - 0.6 control[4] + control[23] <= 0.18
 R85: 0.6 control[4] + control[23] >= -0.18
 R86: control[42] - mode[4] >= -1
 R87: control[42] + mode[4] <= 1
 R88: control[42] - 5 mode[23] >= -5
 R89: 0.6 control[4] - control[23] - mode[23] >= -1.18
 R90: 0.6 control[4] - control[23] + mode[23] <= 0.82
 R91: control[42] + 5 mode[42] <= 5
 R92: 0.6 control[4] + control[23] - mode[42] >= -1.18
 R93: 0.6 control[4] + control[23] + mode[42] <= 0.82
 R94: mode[4] + mode[23] + mode[42] = 1
 R95: - state[5] + state[6] - 0.00261806778230845 control[5] = 0
 R96: - state[25] + state[26] - 7.8542033469253367e-04 state[45]
   - 0.00261806778230845 control[24] = 0
 R97: - state[45] + state[46] + 0.0104166666666667 state[65]
   + 0.0347222222222222 control[24] = 0
 R98: - state[65] + state[66] - 0.05 control[43] = 0
 R99: 0.0523613556461689 control[5] >= -0.0157084066938507
 R100: 0.009375 state[65] + 0.0836113556461689 control[24]
   + 0.045 control[43] >= -0.3
 R101: 0.0523613556461689 control[5] <= 0.2842915933061493
 R102: 0.009375 state[65] + 0.0836113556461689 control[24]
   + 0.045 control[43] <= 0.3
 R103: - 0.6 control[5] + control[24] <= 0.18
 R104: 0.6 control[5] + control[24] >= -0.18
 R105: control[43] - mode[5] >= -1
 R106: control[43] + mode[5] <= 1
 R107: control[43] - 5 mode[24] >= -5
 R108: 0.6 control[5] - control[24] - mode[24] >= -1.18
 R109: 0.6 control[5] - control[24] + mode[24] <= 0.82
 R110: control[43] + 5 mode[43] <= 5
 R111: 0.6 control[5] + control[24] - mode[43] >= -1.18
 R112: 0.6 control[5] + control[24] + mode[43] <= 0.82
 R113: mode[5] + mode[24] + mode[43] = 1
 R114: - state[6] + state[7] - 0.00261806778230845 control[6] = 0
 R115: - state[26] + state[27] - 7.8542033469253367e-04 state[46]
   - 0.00261806778230845 control[25] = 0
 R116: - state[46] + state[47] + 0.0104166666666667 state[66]
   + 0.0347222222222222 control[25] = 0
 R117: - state[66] + state[67] - 0.05 control[44] = 0
 R118: 0.0523613556461689 control[6] >= -0.0157084066938507
 R119: 0.009375 state[66] + 0.0836113556461689 control[25]
   + 0.045 control[44] >= -0.3
 R120: 0.0523613556461689 control[6] <= 0.2842915933061493
 R121: 0.009375 state[66] + 0.0836113556461689 control[25]
   + 0.045 control[44] <= 0.3
 R122: - 0.6 control[6] + control[25] <= 0.18
 R123: 0.6 control[6] + control[25] >= -0.18
 R124: control[44] - mode[6] >= -1
 R125: control[44] + mode[6] <= 1
 R126: control[44] - 5 mode[25] >= -5
 R127: 0.6 control[6] - control[25] - mode[25] >= -1.18
 R128: 0.6 control[6] - control[25] + mode[25] <= 0.82
 R129: control[44] + 5 mode[44] <= 5
 R130: 0.6 control[6] + control[25] - mode[44] >= -1.18
 R131: 0.6 control[6] + control[25] + mode[44] <= 0.82
 R132: mode[6] + mode[25] + mode[44] = 1
 R133: - state[7] + state[8] - 0.00261806778230845 control[7] = 0
 R134: - state[27] + state[28] - 7.8542033469253367e-04 state[47]
   - 0.00261806778230845 control[26] = 0
 R135: - state[47] + state[48] + 0.0104166666666667 state[67]
   + 0.0347222222222222 control[26] = 0
 R136: - state[67] + state[68] - 0.05 control[45] = 0
 R137: 0.0523613556461689 control[7] >= -0.0157084066938507
 R138: 0.009375 state[67] + 0.0836113556461689 control[26]
   + 0.045 control[45] >= -0.3
 R139: 0.0523613556461689 control[7] <= 0.2842915933061493
 R140: 0.009375 state[67] + 0.0836113556461689 control[26]
   + 0.045 control[45] <= 0.3
 R141: - 0.6 control[7] + control[26] <= 0.18
 R142: 0.6 control[7] + control[26] >= -0.18
 R143: control[45] - mode[7] >= -1
 R144: control[45] + mode[7] <= 1
 R145: control[45] - 5 mode[26] >= -5
 R146: 0.6 control[7] - control[26] - mode[26] >= -1.18
 R147: 0.6 control[7] - control[26] + mode[26] <= 0.82
 R148: control[45] + 5 mode[45] <= 5
 R149: 0.6 control[7] + control[26] - mode[45] >= -1.18
 R150: 0.6 control[7] + control[26] + mode[45] <= 0.82
 R151: mode[7] + mode[26] + mode[45] = 1
 R152: - state[8] + state[9] - 0.00261806778230845 control[8] = 0
 R153: - state[28] + state[29] - 7.8542033469253367e-04 state[48]
   - 0.00261806778230845 control[27] = 0
 R154: - state[48] + state[49] + 0.0104166666666667 state[68]
   + 0.0347222222222222 control[27] = 0
 R155: - state[68] + state[69] - 0.05 control[46] = 0
 R156: 0.0523613556461689 control[8] >= -0.0157084066938507
 R157: 0.009375 state[68] + 0.0836113556461689 control[27]
   + 0.045 control[46] >= -0.3
 R158: 0.0523613556461689 control[8] <= 0.2842915933061493
 R159: 0.009375 state[68] + 0.0836113556461689 control[27]
   + 0.045 control[46] <= 0.3
 R160: - 0.6 control[8] + control[27] <= 0.18
 R161: 0.6 control[8] + control[27] >= -0.18
 R162: control[46] - mode[8] >= -1
 R163: control[46] + mode[8] <= 1
 R164: control[46] - 5 mode[27] >= -5
 R165: 0.6 control[8] - control[27] - mode[27] >= -1.18
 R166: 0.6 control[8] - control[27] + mode[27] <= 0.82
 R167: control[46] + 5 mode[46] <= 5
 R168: 0.6 control[8] + control[27] - mode[46] >= -1.18
 R169: 0.6 control[8] + control[27] + mode[46] <= 0.82
 R170: mode[8] + mode[27] + mode[46] = 1
 R171: - state[9] + state[10] - 0.00261806778230845 control[9] = 0
 R172: - state[29] + state[30] - 7.8542033469253367e-04 state[49]
   - 0.00261806778230845 control[28] = 0
 R173: - state[49] + state[50] + 0.0104166666666667 state[69]
   + 0.0347222222222222 control[28] = 0
 R174: - state[69] + state[70] - 0.05 control[47] = 0
 R175: 0.0523613556461689 control[9] >= -0.0157084066938507
 R176: 0.009375 state[69] + 0.0836113556461689 control[28]
   + 0.045 control[47] >= -0.3
 R177: 0.0523613556461689 control[9] <= 0.2842915933061493
 R178: 0.009375 state[69] + 0.0836113556461689 control[28]
   + 0.045 control[47] <= 0.3
 R179: - 0.6 control[9] + control[28] <= 0.18
 R180: 0.6 control[9] + control[28] >= -0.18
 R181: control[47] - mode[9] >= -1
 R182: control[47] + mode[9] <= 1
 R183: control[47] - 5 mode[28] >= -5
 R184: 0.6 control[9] - control[28] - mode[28] >= -1.18
 R185: 0.6 control[9] - control[28] + mode[28] <= 0.82
 R186: control[47] + 5 mode[47] <= 5
 R187: 0.6 control[9] + control[28] - mode[47] >= -1.18
 R188: 0.6 control[9] + control[28] + mode[47] <= 0.82
 R189: mode[9] + mode[28] + mode[47] = 1
 R190: - state[10] + state[11] - 0.00261806778230845 control[10] = 0
 R191: - state[30] + state[31] - 7.8542033469253367e-04 state[50]
   - 0.00261806778230845 control[29] = 0
 R192: - state[50] + state[51] + 0.0104166666666667 state[70]
   + 0.0347222222222222 control[29] = 0
 R193: - state[70] + state[71] - 0.05 control[48] = 0
 R194: 0.0523613556461689 control[10] >= -0.0157084066938507
 R195: 0.009375 state[70] + 0.0836113556461689 control[29]
   + 0.045 control[48] >= -0.3
 R196: 0.0523613556461689 control[10] <= 0.2842915933061493
 R197: 0.009375 state[70] + 0.0836113556461689 control[29]
   + 0.045 control[48] <= 0.3
 R198: - 0.6 control[10] + control[29] <= 0.18
 R199: 0.6 control[10] + control[29] >= -0.18
 R200: control[48] - mode[10] >= -1
 R201: control[48] + mode[10] <= 1
 R202: control[48] - 5 mode[29] >= -5
 R203: 0.6 control[10] - control[29] - mode[29] >= -1.18
 R204: 0.6 control[10] - control[29] + mode[29] <= 0.82
 R205: control[48] + 5 mode[48] <= 5
 R206: 0.6 control[10] + control[29] - mode[48] >= -1.18
 R207: 0.6 control[10] + control[29] + mode[48] <= 0.82
 R208: mode[10] + mode[29] + mode[48] = 1
 R209: - state[11] + state[12] - 0.00261806778230845 control[11] = 0
 R210: - state[31] + state[32] - 7.8542033469253367e-04 state[51]
   - 0.00261806778230845 control[30] = 0
 R211: - state[51] + state[52] + 0.0104166666666667 state[71]
   + 0.0347222222222222 control[30] = 0
 R212: - state[71] + state[72] - 0.05 control[49] = 0
 R213: 0.0523613556461689 control[11] >= -0.0157084066938507
 R214: 0.009375 state[71] + 0.0836113556461689 control[30]
   + 0.045 control[49] >= -0.3
 R215: 0.0523613556461689 control[11] <= 0.2842915933061493
 R216: 0.009375 state[71] + 0.0836113556461689 control[30]
   + 0.045 control[49] <= 0.3
 R217: - 0.6 control[11] + control[30] <= 0.18
 R218: 0.6 control[11] + control[30] >= -0.18
 R219: control[49] - mode[11] >= -1
 R220: control[49] + mode[11] <= 1
 R221: control[49] - 5 mode[30] >= -5
 R222: 0.6 control[11] - control[30] - mode[30] >= -1.18
 R223: 0.6 control[11] - control[30] + mode[30] <= 0.82
 R224: control[49] + 5 mode[49] <= 5
 R225: 0.6 control[11] + control[30] - mode[49] >= -1.18
 R226: 0.6 control[11] + control[30] + mode[49] <= 0.82
 R227: mode[11] + mode[30] + mode[49] = 1
 R228: - state[12] + state[13] - 0.00261806778230845 control[12] = 0
 R229: - state[32] + state[33] - 7.8542033469253367e-04 state[52]
   - 0.00261806778230845 control[31] = 0
 R230: - state[52] + state[53] + 0.0104166666666667 state[72]
   + 0.0347222222222222 control[31] = 0
 R231: - state[72] + state[73] - 0.05 control[50] = 0
 R232: 0.0523613556461689 control[12] >= -0.0157084066938507
 R233: 0.009375 state[72] + 0.0836113556461689 control[31]
   + 0.045 control[50] >= -0.3
 R234: 0.0523613556461689 control[12] <= 0.2842915933061493
 R235: 0.009375 state[72] + 0.0836113556461689 control[31]
   + 0.045 control[50] <= 0.3
 R236: - 0.6 control[12] + control[31] <= 0.18
 R237: 0.6 control[12] + control[31] >= -0.18
 R238: control[50] - mode[12] >= -1
 R239: control[50] + mode[12] <= 1
 R240: control[50] - 5 mode[31] >= -5
 R241: 0.6 control[12] - control[31] - mode[31] >= -1.18
 R242: 0.6 control[12] - control[31] + mode[31] <= 0.82
 R243: control[50] + 5 mode[50] <= 5
 R244: 0.6 control[12] + control[31] - mode[50] >= -1.18
 R245: 0.6 control[12] + control[31] + mode[50] <= 0.82
 R246: mode[12] + mode[31] + mode[50] = 1
 R247: - state[13] + state[14] - 0.00261806778230845 control[13] = 0
 R248: - state[33] + state[34] - 7.8542033469253367e-04 state[53]
   - 0.00261806778230845 control[32] = 0
 R249: - state[53] + state[54] + 0.0104166666666667 state[73]
   + 0.0347222222222222 control[32] = 0
 R250: - state[73] + state[74] - 0.05 control[51] = 0
 R251: 0.0523613556461689 control[13] >= -0.0157084066938507
 R252: 0.009375 state[73] + 0.0836113556461689 control[32]
   + 0.045 control[51] >= -0.3
 R253: 0.0523613556461689 control[13] <= 0.2842915933061493
 R254: 0.009375 state[73] + 0.0836113556461689 control[32]
   + 0.045 control[51] <= 0.3
 R255: - 0.6 control[13] + control[32] <= 0.18
 R256: 0.6 control[13] + control[32] >= -0.18
 R257: control[51] - mode[13] >= -1
 R258: control[51] + mode[13] <= 1
 R259: control[51] - 5 mode[32] >= -5
 R260: 0.6 control[13] - control[32] - mode[32] >= -1.18
 R261: 0.6 control[13] - control[32] + mode[32] <= 0.82
 R262: control[51] + 5 mode[51] <= 5
 R263: 0.6 control[13] + control[32] - mode[51] >= -1.18
 R264: 0.6 control[13] + control[32] + mode[51] <= 0.82
 R265: mode[13] + mode[32] + mode[51] = 1
 R266: - state[14] + state[15] - 0.00261806778230845 control[14] = 0
 R267: - state[34] + state[35] - 7.8542033469253367e-04 state[54]
   - 0.00261806778230845 control[33] = 0
 R268: - state[54] + state[55] + 0.0104166666666667 state[74]
   + 0.0347222222222222 control[33] = 0
 R269: - state[74] + state[75] - 0.05 control[52] = 0
 R270: 0.0523613556461689 control[14] >= -0.0157084066938507
 R271: 0.009375 state[74] + 0.0836113556461689 control[33]
   + 0.045 control[52] >= -0.3
 R272: 0.0523613556461689 control[14] <= 0.2842915933061493
 R273: 0.009375 state[74] + 0.0836113556461689 control[33]
   + 0.045 control[52] <= 0.3
 R274: - 0.6 control[14] + control[33] <= 0.18
 R275: 0.6 control[14] + control[33] >= -0.18
 R276: control[52] - mode[14] >= -1
 R277: control[52] + mode[14] <= 1
 R278: control[52] - 5 mode[33] >= -5
 R279: 0.6 control[14] - control[33] - mode[33] >= -1.18
 R280: 0.6 control[14] - control[33] + mode[33] <= 0.82
 R281: control[52] + 5 mode[52] <= 5
 R282: 0.6 control[14] + control[33] - mode[52] >= -1.18
 R283: 0.6 control[14] + control[33] + mode[52] <= 0.82
 R284: mode[14] + mode[33] + mode[52] = 1
 R285: - state[15] + state[16] - 0.00261806778230845 control[15] = 0
 R286: - state[35] + state[36] - 7.8542033469253367e-04 state[55]
   - 0.00261806778230845 control[34] = 0
 R287: - state[55] + state[56] + 0.0104166666666667 state[75]
   + 0.0347222222222222 control[34] = 0
 R288: - state[75] + state[76] - 0.05 control[53] = 0
 R289: 0.0523613556461689 control[15] >= -0.0157084066938507
 R290: 0.009375 state[75] + 0.0836113556461689 control[34]
   + 0.045 control[53] >= -0.3
 R291: 0.0523613556461689 control[15] <= 0.2842915933061493
 R292: 0.009375 state[75] + 0.0836113556461689 control[34]
   + 0.045 control[53] <= 0.3
 R293: - 0.6 control[15] + control[34] <= 0.18
 R294: 0.6 control[15] + control[34] >= -0.18
 R295: control[53] - mode[15] >= -1
 R296: control[53] + mode[15] <= 1
 R297: control[53] - 5 mode[34] >= -5
 R298: 0.6 control[15] - control[34] - mode[34] >= -1.18
 R299: 0.6 control[15] - control[34] + mode[34] <= 0.82
 R300: control[53] + 5 mode[53] <= 5
 R301: 0.6 control[15] + control[34] - mode[53] >= -1.18
 R302: 0.6 control[15] + control[34] + mode[53] <= 0.82
 R303: mode[15] + mode[34] + mode[53] = 1
 R304: - state[16] + state[17] - 0.00261806778230845 control[16] = 0
 R305: - state[36] + state[37] - 7.8542033469253367e-04 state[56]
   - 0.00261806778230845 control[35] = 0
 R306: - state[56] + state[57] + 0.0104166666666667 state[76]
   + 0.0347222222222222 control[35] = 0
 R307: - state[76] + state[77] - 0.05 control[54] = 0
 R308: 0.0523613556461689 control[16] >= -0.0157084066938507
 R309: 0.009375 state[76] + 0.0836113556461689 control[35]
   + 0.045 control[54] >= -0.3
 R310: 0.0523613556461689 control[16] <= 0.2842915933061493
 R311: 0.009375 state[76] + 0.0836113556461689 control[35]
   + 0.045 control[54] <= 0.3
 R312: - 0.6 control[16] + control[35] <= 0.18
 R313: 0.6 control[16] + control[35] >= -0.18
 R314: control[54] - mode[16] >= -1
 R315: control[54] + mode[16] <= 1
 R316: control[54] - 5 mode[35] >= -5
 R317: 0.6 control[16] - control[35] - mode[35] >= -1.18
 R318: 0.6 control[16] - control[35] + mode[35] <= 0.82
 R319: control[54] + 5 mode[54] <= 5
 R320: 0.6 control[16] + control[35] - mode[54] >= -1.18
 R321: 0.6 control[16] + control[35] + mode[54] <= 0.82
 R322: mode[16] + mode[35] + mode[54] = 1
 R323: - state[17] + state[18] - 0.00261806778230845 control[17] = 0
 R324: - state[37] + state[38] - 7.8542033469253367e-04 state[57]
   - 0.00261806778230845 control[36] = 0
 R325: - state[57] + state[58] + 0.0104166666666667 state[77]
   + 0.0347222222222222 control[36] = 0
 R326: - state[77] + state[78] - 0.05 control[55] = 0
 R327: 0.0523613556461689 control[17] >= -0.0157084066938507
 R328: 0.009375 state[77] + 0.0836113556461689 control[36]
   + 0.045 control[55] >= -0.3
 R329: 0.0523613556461689 control[17] <= 0.2842915933061493
 R330: 0.009375 state[77] + 0.0836113556461689 control[36]
   + 0.045 control[55] <= 0.3
 R331: - 0.6 control[17] + control[36] <= 0.18
 R332: 0.6 control[17] + control[36] >= -0.18
 R333: control[55] - mode[17] >= -1
 R334: control[55] + mode[17] <= 1
 R335: control[55] - 5 mode[36] >= -5
 R336: 0.6 control[17] - control[36] - mode[36] >= -1.18
 R337: 0.6 control[17] - control[36] + mode[36] <= 0.82
 R338: control[55] + 5 mode[55] <= 5
 R339: 0.6 control[17] + control[36] - mode[55] >= -1.18
 R340: 0.6 control[17] + control[36] + mode[55] <= 0.82
 R341: mode[17] + mode[36] + mode[55] = 1
 R342: - state[18] + state[19] - 0.00261806778230845 control[18] = 0
 R343: - state[38] + state[39] - 7.8542033469253367e-04 state[58]
   - 0.00261806778230845 control[37] = 0
 R344: - state[58] + state[59] + 0.0104166666666667 state[78]
   + 0.0347222222222222 control[37] = 0
 R345: - state[78] + state[79] - 0.05 control[56] = 0
 R346: 0.0523613556461689 control[18] >= -0.0157084066938507
 R347: 0.009375 state[78] + 0.0836113556461689 control[37]
   + 0.045 control[56] >= -0.3
 R348: 0.0523613556461689 control[18] <= 0.2842915933061493
 R349: 0.009375 state[78] + 0.0836113556461689 control[37]
   + 0.045 control[56] <= 0.3
 R350: - 0.6 control[18] + control[37] <= 0.18
 R351: 0.6 control[18] + control[37] >= -0.18
 R352: control[56] - mode[18] >= -1
 R353: control[56] + mode[18] <= 1
 R354: control[56] - 5 mode[37] >= -5
 R355: 0.6 control[18] - control[37] - mode[37] >= -1.18
 R356: 0.6 control[18] - control[37] + mode[37] <= 0.82
 R357: control[56] + 5 mode[56] <= 5
 R358: 0.6 control[18] + control[37] - mode[56] >= -1.18
 R359: 0.6 control[18] + control[37] + mode[56] <= 0.82
 R360: mode[18] + mode[37] + mode[56] = 1
Bounds
 state[0] = 0
 state[1] free
 state[2] free
 state[3] free
 state[4] free
 state[5] free
 state[6] free
 state[7] free
 state[8] free
 state[9] free
 state[10] free
 state[11] free
 state[12] free
 state[13] free
 state[14] free
 state[15] free
 state[16] free
 state[17] free
 state[18] free
 state[19] free
 state[20] = 0.04
 state[21] free
 state[22] free
 state[23] free
 state[24] free
 state[25] free
 state[26] free
 state[27] free
 state[28] free
 state[29] free
 state[30] free
 state[31] free
 state[32] free
 state[33] free
 state[34] free
 state[35] free
 state[36] free
 state[37] free
 state[38] free
 state[39] free
 state[40] = 0.06
 state[41] free
 state[42] free
 state[43] free
 state[44] free
 state[45] free
 state[46] free
 state[47] free
 state[48] free
 state[49] free
 state[50] free
 state[51] free
 state[52] free
 state[53] free
 state[54] free
 state[55] free
 state[56] free
 state[57] free
 state[58] free
 state[59] free
 state[60] = -0.5
 -0.7 <= state[61] <= 0.7
 -0.7 <= state[62] <= 0.7
 -0.7 <= state[63] <= 0.7
 -0.7 <= state[64] <= 0.7
 -0.7 <= state[65] <= 0.7
 -0.7 <= state[66] <= 0.7
 -0.7 <= state[67] <= 0.7
 -0.7 <= state[68] <= 0.7
 -0.7 <= state[69] <= 0.7
 -0.7 <= state[70] <= 0.7
 -0.7 <= state[71] <= 0.7
 -0.7 <= state[72] <= 0.7
 -0.7 <= state[73] <= 0.7
 -0.7 <= state[74] <= 0.7
 -0.7 <= state[75] <= 0.7
 -0.7 <= state[76] <= 0.7
 -0.7 <= state[77] <= 0.7
 -0.7 <= state[78] <= 0.7
 -0.7 <= state[79] <= 0.7
 control[0] >= -0.3
 control[1] >= -0.3
 control[2] >= -0.3
 control[3] >= -0.3
 control[4] >= -0.3
 control[5] >= -0.3
 control[6] >= -0.3
 control[7] >= -0.3
 control[8] >= -0.3
 control[9] >= -0.3
 control[10] >= -0.3
 control[11] >= -0.3
 control[12] >= -0.3
 control[13] >= -0.3
 control[14] >= -0.3
 control[15] >= -0.3
 control[16] >= -0.3
 control[17] >= -0.3
 control[18] >= -0.3
 control[19] free
 control[20] free
 control[21] free
 control[22] free
 control[23] free
 control[24] free
 control[25] free
 control[26] free
 control[27] free
 control[28] free
 control[29] free
 control[30] free
 control[31] free
 control[32] free
 control[33] free
 control[34] free
 control[35] free
 control[36] free
 control[37] free
 control[38] free
 control[39] free
 control[40] free
 control[41] free
 control[42] free
 control[43] free
 control[44] free
 control[45] free
 control[46] free
 control[47] free
 control[48] free
 control[49] free
 control[50] free
 control[51] free
 control[52] free
 control[53] free
 control[54] free
 control[55] free
 control[56] free
Binaries
 mode[0] mode[1] mode[2] mode[3] mode[4] mode[5] mode[6] mode[7] mode[8]
 mode[9] mode[10] mode[11] mode[12] mode[13] mode[14] mode[15] mode[16]
 mode[17] mode[18] mode[19] mode[20] mode[21] mode[22] mode[23] mode[24]
 mode[25] mode[26] mode[27] mode[28] mode[29] mode[30] mode[31] mode[32]
 mode[33] mode[34] mode[35] mode[36] mode[37] mode[38] mode[39] mode[40]
 mode[41] mode[42] mode[43] mode[44] mode[45] mode[46] mode[47] mode[48]
 mode[49] mode[50] mode[51] mode[52] mode[53] mode[54] mode[55] mode[56]
End
