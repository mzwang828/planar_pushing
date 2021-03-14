\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
 [ 60 state[0] ^2 + 60 state[20] ^2 + 2 state[40] ^2 + 60 state[1] ^2
   + 60 state[21] ^2 + 2 state[41] ^2 + 60 state[2] ^2 + 60 state[22] ^2
   + 2 state[42] ^2 + 60 state[3] ^2 + 60 state[23] ^2 + 2 state[43] ^2
   + 60 state[4] ^2 + 60 state[24] ^2 + 2 state[44] ^2 + 60 state[5] ^2
   + 60 state[25] ^2 + 2 state[45] ^2 + 60 state[6] ^2 + 60 state[26] ^2
   + 2 state[46] ^2 + 60 state[7] ^2 + 60 state[27] ^2 + 2 state[47] ^2
   + 60 state[8] ^2 + 60 state[28] ^2 + 2 state[48] ^2 + 60 state[9] ^2
   + 60 state[29] ^2 + 2 state[49] ^2 + 60 state[10] ^2 + 60 state[30] ^2
   + 2 state[50] ^2 + 60 state[11] ^2 + 60 state[31] ^2 + 2 state[51] ^2
   + 60 state[12] ^2 + 60 state[32] ^2 + 2 state[52] ^2 + 60 state[13] ^2
   + 60 state[33] ^2 + 2 state[53] ^2 + 60 state[14] ^2 + 60 state[34] ^2
   + 2 state[54] ^2 + 60 state[15] ^2 + 60 state[35] ^2 + 2 state[55] ^2
   + 60 state[16] ^2 + 60 state[36] ^2 + 2 state[56] ^2 + 60 state[17] ^2
   + 60 state[37] ^2 + 2 state[57] ^2 + 60 state[18] ^2 + 60 state[38] ^2
   + 2 state[58] ^2 + 12000 state[19] ^2 + 12000 state[39] ^2
   + 400 state[59] ^2 + control[0] ^2 + control[19] ^2
   + 0.01 control[38] ^2 + control[1] ^2 + control[20] ^2
   + 0.01 control[39] ^2 + control[2] ^2 + control[21] ^2
   + 0.01 control[40] ^2 + control[3] ^2 + control[22] ^2
   + 0.01 control[41] ^2 + control[4] ^2 + control[23] ^2
   + 0.01 control[42] ^2 + control[5] ^2 + control[24] ^2
   + 0.01 control[43] ^2 + control[6] ^2 + control[25] ^2
   + 0.01 control[44] ^2 + control[7] ^2 + control[26] ^2
   + 0.01 control[45] ^2 + control[8] ^2 + control[27] ^2
   + 0.01 control[46] ^2 + control[9] ^2 + control[28] ^2
   + 0.01 control[47] ^2 + control[10] ^2 + control[29] ^2
   + 0.01 control[48] ^2 + control[11] ^2 + control[30] ^2
   + 0.01 control[49] ^2 + control[12] ^2 + control[31] ^2
   + 0.01 control[50] ^2 + control[13] ^2 + control[32] ^2
   + 0.01 control[51] ^2 + control[14] ^2 + control[33] ^2
   + 0.01 control[52] ^2 + control[15] ^2 + control[34] ^2
   + 0.01 control[53] ^2 + control[16] ^2 + control[35] ^2
   + 0.01 control[54] ^2 + control[17] ^2 + control[36] ^2
   + 0.01 control[55] ^2 + control[18] ^2 + control[37] ^2
   + 0.01 control[56] ^2 ] / 2 
Subject To
 R0: state[0] - state[1] + 0.1538782696540474 control[0] = 0
 R1: state[20] + 0.0469328722444845 state[40] - state[21]
   + 0.1538782696540474 control[19] = 0
 R2: state[40] - 0.2118055555555556 state[60] - state[41]
   - 0.6944444444444444 control[19] = 0
 R3: state[60] - state[61] + control[38] = 0
 R4: 0.1538782696540474 control[0] >= -0.0469328722444845
 R5: 0.1538782696540474 control[0] <= 0.2530671277555155
 R6: 0.00953125 state[60] + 0.1851282696540474 control[19]
   + 0.045 control[38] >= -0.3
 R7: 0.00953125 state[60] + 0.1851282696540474 control[19]
   + 0.045 control[38] <= 0.3
 R8: - 0.6 control[0] + control[19] <= 0.183
 R9: 0.6 control[0] + control[19] >= -0.183
 R10: control[38] - mode[0] >= -1
 R11: control[38] + mode[0] <= 1
 R12: control[38] - 5 mode[19] >= -5
 R13: 0.6 control[0] - control[19] - mode[19] >= -1.183
 R14: 0.6 control[0] - control[19] + mode[19] <= 0.817
 R15: control[38] + 5 mode[38] <= 5
 R16: 0.6 control[0] + control[19] - mode[38] >= -1.183
 R17: 0.6 control[0] + control[19] + mode[38] <= 0.817
 R18: mode[0] + mode[19] + mode[38] = 1
 R19: state[1] - state[2] + 0.1538782696540474 control[1] = 0
 R20: state[21] + 0.0469328722444845 state[41] - state[22]
   + 0.1538782696540474 control[20] = 0
 R21: state[41] - 0.2118055555555556 state[61] - state[42]
   - 0.6944444444444444 control[20] = 0
 R22: state[61] - state[62] + control[39] = 0
 R23: 0.1538782696540474 control[1] >= -0.0469328722444845
 R24: 0.1538782696540474 control[1] <= 0.2530671277555155
 R25: 0.00953125 state[61] + 0.1851282696540474 control[20]
   + 0.045 control[39] >= -0.3
 R26: 0.00953125 state[61] + 0.1851282696540474 control[20]
   + 0.045 control[39] <= 0.3
 R27: - 0.6 control[1] + control[20] <= 0.183
 R28: 0.6 control[1] + control[20] >= -0.183
 R29: control[39] - mode[1] >= -1
 R30: control[39] + mode[1] <= 1
 R31: control[39] - 5 mode[20] >= -5
 R32: 0.6 control[1] - control[20] - mode[20] >= -1.183
 R33: 0.6 control[1] - control[20] + mode[20] <= 0.817
 R34: control[39] + 5 mode[39] <= 5
 R35: 0.6 control[1] + control[20] - mode[39] >= -1.183
 R36: 0.6 control[1] + control[20] + mode[39] <= 0.817
 R37: mode[1] + mode[20] + mode[39] = 1
 R38: state[2] - state[3] + 0.1538782696540474 control[2] = 0
 R39: state[22] + 0.0469328722444845 state[42] - state[23]
   + 0.1538782696540474 control[21] = 0
 R40: state[42] - 0.2118055555555556 state[62] - state[43]
   - 0.6944444444444444 control[21] = 0
 R41: state[62] - state[63] + control[40] = 0
 R42: 0.1538782696540474 control[2] >= -0.0469328722444845
 R43: 0.1538782696540474 control[2] <= 0.2530671277555155
 R44: 0.00953125 state[62] + 0.1851282696540474 control[21]
   + 0.045 control[40] >= -0.3
 R45: 0.00953125 state[62] + 0.1851282696540474 control[21]
   + 0.045 control[40] <= 0.3
 R46: - 0.6 control[2] + control[21] <= 0.183
 R47: 0.6 control[2] + control[21] >= -0.183
 R48: control[40] - mode[2] >= -1
 R49: control[40] + mode[2] <= 1
 R50: control[40] - 5 mode[21] >= -5
 R51: 0.6 control[2] - control[21] - mode[21] >= -1.183
 R52: 0.6 control[2] - control[21] + mode[21] <= 0.817
 R53: control[40] + 5 mode[40] <= 5
 R54: 0.6 control[2] + control[21] - mode[40] >= -1.183
 R55: 0.6 control[2] + control[21] + mode[40] <= 0.817
 R56: mode[2] + mode[21] + mode[40] = 1
 R57: state[3] - state[4] + 0.1538782696540474 control[3] = 0
 R58: state[23] + 0.0469328722444845 state[43] - state[24]
   + 0.1538782696540474 control[22] = 0
 R59: state[43] - 0.2118055555555556 state[63] - state[44]
   - 0.6944444444444444 control[22] = 0
 R60: state[63] - state[64] + control[41] = 0
 R61: 0.1538782696540474 control[3] >= -0.0469328722444845
 R62: 0.1538782696540474 control[3] <= 0.2530671277555155
 R63: 0.00953125 state[63] + 0.1851282696540474 control[22]
   + 0.045 control[41] >= -0.3
 R64: 0.00953125 state[63] + 0.1851282696540474 control[22]
   + 0.045 control[41] <= 0.3
 R65: - 0.6 control[3] + control[22] <= 0.183
 R66: 0.6 control[3] + control[22] >= -0.183
 R67: control[41] - mode[3] >= -1
 R68: control[41] + mode[3] <= 1
 R69: control[41] - 5 mode[22] >= -5
 R70: 0.6 control[3] - control[22] - mode[22] >= -1.183
 R71: 0.6 control[3] - control[22] + mode[22] <= 0.817
 R72: control[41] + 5 mode[41] <= 5
 R73: 0.6 control[3] + control[22] - mode[41] >= -1.183
 R74: 0.6 control[3] + control[22] + mode[41] <= 0.817
 R75: mode[3] + mode[22] + mode[41] = 1
 R76: state[4] - state[5] + 0.1538782696540474 control[4] = 0
 R77: state[24] + 0.0469328722444845 state[44] - state[25]
   + 0.1538782696540474 control[23] = 0
 R78: state[44] - 0.2118055555555556 state[64] - state[45]
   - 0.6944444444444444 control[23] = 0
 R79: state[64] - state[65] + control[42] = 0
 R80: 0.1538782696540474 control[4] >= -0.0469328722444845
 R81: 0.1538782696540474 control[4] <= 0.2530671277555155
 R82: 0.00953125 state[64] + 0.1851282696540474 control[23]
   + 0.045 control[42] >= -0.3
 R83: 0.00953125 state[64] + 0.1851282696540474 control[23]
   + 0.045 control[42] <= 0.3
 R84: - 0.6 control[4] + control[23] <= 0.183
 R85: 0.6 control[4] + control[23] >= -0.183
 R86: control[42] - mode[4] >= -1
 R87: control[42] + mode[4] <= 1
 R88: control[42] - 5 mode[23] >= -5
 R89: 0.6 control[4] - control[23] - mode[23] >= -1.183
 R90: 0.6 control[4] - control[23] + mode[23] <= 0.817
 R91: control[42] + 5 mode[42] <= 5
 R92: 0.6 control[4] + control[23] - mode[42] >= -1.183
 R93: 0.6 control[4] + control[23] + mode[42] <= 0.817
 R94: mode[4] + mode[23] + mode[42] = 1
 R95: state[5] - state[6] + 0.1538782696540474 control[5] = 0
 R96: state[25] + 0.0469328722444845 state[45] - state[26]
   + 0.1538782696540474 control[24] = 0
 R97: state[45] - 0.2118055555555556 state[65] - state[46]
   - 0.6944444444444444 control[24] = 0
 R98: state[65] - state[66] + control[43] = 0
 R99: 0.1538782696540474 control[5] >= -0.0469328722444845
 R100: 0.1538782696540474 control[5] <= 0.2530671277555155
 R101: 0.00953125 state[65] + 0.1851282696540474 control[24]
   + 0.045 control[43] >= -0.3
 R102: 0.00953125 state[65] + 0.1851282696540474 control[24]
   + 0.045 control[43] <= 0.3
 R103: - 0.6 control[5] + control[24] <= 0.183
 R104: 0.6 control[5] + control[24] >= -0.183
 R105: control[43] - mode[5] >= -1
 R106: control[43] + mode[5] <= 1
 R107: control[43] - 5 mode[24] >= -5
 R108: 0.6 control[5] - control[24] - mode[24] >= -1.183
 R109: 0.6 control[5] - control[24] + mode[24] <= 0.817
 R110: control[43] + 5 mode[43] <= 5
 R111: 0.6 control[5] + control[24] - mode[43] >= -1.183
 R112: 0.6 control[5] + control[24] + mode[43] <= 0.817
 R113: mode[5] + mode[24] + mode[43] = 1
 R114: state[6] - state[7] + 0.1538782696540474 control[6] = 0
 R115: state[26] + 0.0469328722444845 state[46] - state[27]
   + 0.1538782696540474 control[25] = 0
 R116: state[46] - 0.2118055555555556 state[66] - state[47]
   - 0.6944444444444444 control[25] = 0
 R117: state[66] - state[67] + control[44] = 0
 R118: 0.1538782696540474 control[6] >= -0.0469328722444845
 R119: 0.1538782696540474 control[6] <= 0.2530671277555155
 R120: 0.00953125 state[66] + 0.1851282696540474 control[25]
   + 0.045 control[44] >= -0.3
 R121: 0.00953125 state[66] + 0.1851282696540474 control[25]
   + 0.045 control[44] <= 0.3
 R122: - 0.6 control[6] + control[25] <= 0.183
 R123: 0.6 control[6] + control[25] >= -0.183
 R124: control[44] - mode[6] >= -1
 R125: control[44] + mode[6] <= 1
 R126: control[44] - 5 mode[25] >= -5
 R127: 0.6 control[6] - control[25] - mode[25] >= -1.183
 R128: 0.6 control[6] - control[25] + mode[25] <= 0.817
 R129: control[44] + 5 mode[44] <= 5
 R130: 0.6 control[6] + control[25] - mode[44] >= -1.183
 R131: 0.6 control[6] + control[25] + mode[44] <= 0.817
 R132: mode[6] + mode[25] + mode[44] = 1
 R133: state[7] - state[8] + 0.1538782696540474 control[7] = 0
 R134: state[27] + 0.0469328722444845 state[47] - state[28]
   + 0.1538782696540474 control[26] = 0
 R135: state[47] - 0.2118055555555556 state[67] - state[48]
   - 0.6944444444444444 control[26] = 0
 R136: state[67] - state[68] + control[45] = 0
 R137: 0.1538782696540474 control[7] >= -0.0469328722444845
 R138: 0.1538782696540474 control[7] <= 0.2530671277555155
 R139: 0.00953125 state[67] + 0.1851282696540474 control[26]
   + 0.045 control[45] >= -0.3
 R140: 0.00953125 state[67] + 0.1851282696540474 control[26]
   + 0.045 control[45] <= 0.3
 R141: - 0.6 control[7] + control[26] <= 0.183
 R142: 0.6 control[7] + control[26] >= -0.183
 R143: control[45] - mode[7] >= -1
 R144: control[45] + mode[7] <= 1
 R145: control[45] - 5 mode[26] >= -5
 R146: 0.6 control[7] - control[26] - mode[26] >= -1.183
 R147: 0.6 control[7] - control[26] + mode[26] <= 0.817
 R148: control[45] + 5 mode[45] <= 5
 R149: 0.6 control[7] + control[26] - mode[45] >= -1.183
 R150: 0.6 control[7] + control[26] + mode[45] <= 0.817
 R151: mode[7] + mode[26] + mode[45] = 1
 R152: state[8] - state[9] + 0.1538782696540474 control[8] = 0
 R153: state[28] + 0.0469328722444845 state[48] - state[29]
   + 0.1538782696540474 control[27] = 0
 R154: state[48] - 0.2118055555555556 state[68] - state[49]
   - 0.6944444444444444 control[27] = 0
 R155: state[68] - state[69] + control[46] = 0
 R156: 0.1538782696540474 control[8] >= -0.0469328722444845
 R157: 0.1538782696540474 control[8] <= 0.2530671277555155
 R158: 0.00953125 state[68] + 0.1851282696540474 control[27]
   + 0.045 control[46] >= -0.3
 R159: 0.00953125 state[68] + 0.1851282696540474 control[27]
   + 0.045 control[46] <= 0.3
 R160: - 0.6 control[8] + control[27] <= 0.183
 R161: 0.6 control[8] + control[27] >= -0.183
 R162: control[46] - mode[8] >= -1
 R163: control[46] + mode[8] <= 1
 R164: control[46] - 5 mode[27] >= -5
 R165: 0.6 control[8] - control[27] - mode[27] >= -1.183
 R166: 0.6 control[8] - control[27] + mode[27] <= 0.817
 R167: control[46] + 5 mode[46] <= 5
 R168: 0.6 control[8] + control[27] - mode[46] >= -1.183
 R169: 0.6 control[8] + control[27] + mode[46] <= 0.817
 R170: mode[8] + mode[27] + mode[46] = 1
 R171: state[9] - state[10] + 0.1538782696540474 control[9] = 0
 R172: state[29] + 0.0469328722444845 state[49] - state[30]
   + 0.1538782696540474 control[28] = 0
 R173: state[49] - 0.2118055555555556 state[69] - state[50]
   - 0.6944444444444444 control[28] = 0
 R174: state[69] - state[70] + control[47] = 0
 R175: 0.1538782696540474 control[9] >= -0.0469328722444845
 R176: 0.1538782696540474 control[9] <= 0.2530671277555155
 R177: 0.00953125 state[69] + 0.1851282696540474 control[28]
   + 0.045 control[47] >= -0.3
 R178: 0.00953125 state[69] + 0.1851282696540474 control[28]
   + 0.045 control[47] <= 0.3
 R179: - 0.6 control[9] + control[28] <= 0.183
 R180: 0.6 control[9] + control[28] >= -0.183
 R181: control[47] - mode[9] >= -1
 R182: control[47] + mode[9] <= 1
 R183: control[47] - 5 mode[28] >= -5
 R184: 0.6 control[9] - control[28] - mode[28] >= -1.183
 R185: 0.6 control[9] - control[28] + mode[28] <= 0.817
 R186: control[47] + 5 mode[47] <= 5
 R187: 0.6 control[9] + control[28] - mode[47] >= -1.183
 R188: 0.6 control[9] + control[28] + mode[47] <= 0.817
 R189: mode[9] + mode[28] + mode[47] = 1
 R190: state[10] - state[11] + 0.1538782696540474 control[10] = 0
 R191: state[30] + 0.0469328722444845 state[50] - state[31]
   + 0.1538782696540474 control[29] = 0
 R192: state[50] - 0.2118055555555556 state[70] - state[51]
   - 0.6944444444444444 control[29] = 0
 R193: state[70] - state[71] + control[48] = 0
 R194: 0.1538782696540474 control[10] >= -0.0469328722444845
 R195: 0.1538782696540474 control[10] <= 0.2530671277555155
 R196: 0.00953125 state[70] + 0.1851282696540474 control[29]
   + 0.045 control[48] >= -0.3
 R197: 0.00953125 state[70] + 0.1851282696540474 control[29]
   + 0.045 control[48] <= 0.3
 R198: - 0.6 control[10] + control[29] <= 0.183
 R199: 0.6 control[10] + control[29] >= -0.183
 R200: control[48] - mode[10] >= -1
 R201: control[48] + mode[10] <= 1
 R202: control[48] - 5 mode[29] >= -5
 R203: 0.6 control[10] - control[29] - mode[29] >= -1.183
 R204: 0.6 control[10] - control[29] + mode[29] <= 0.817
 R205: control[48] + 5 mode[48] <= 5
 R206: 0.6 control[10] + control[29] - mode[48] >= -1.183
 R207: 0.6 control[10] + control[29] + mode[48] <= 0.817
 R208: mode[10] + mode[29] + mode[48] = 1
 R209: state[11] - state[12] + 0.1538782696540474 control[11] = 0
 R210: state[31] + 0.0469328722444845 state[51] - state[32]
   + 0.1538782696540474 control[30] = 0
 R211: state[51] - 0.2118055555555556 state[71] - state[52]
   - 0.6944444444444444 control[30] = 0
 R212: state[71] - state[72] + control[49] = 0
 R213: 0.1538782696540474 control[11] >= -0.0469328722444845
 R214: 0.1538782696540474 control[11] <= 0.2530671277555155
 R215: 0.00953125 state[71] + 0.1851282696540474 control[30]
   + 0.045 control[49] >= -0.3
 R216: 0.00953125 state[71] + 0.1851282696540474 control[30]
   + 0.045 control[49] <= 0.3
 R217: - 0.6 control[11] + control[30] <= 0.183
 R218: 0.6 control[11] + control[30] >= -0.183
 R219: control[49] - mode[11] >= -1
 R220: control[49] + mode[11] <= 1
 R221: control[49] - 5 mode[30] >= -5
 R222: 0.6 control[11] - control[30] - mode[30] >= -1.183
 R223: 0.6 control[11] - control[30] + mode[30] <= 0.817
 R224: control[49] + 5 mode[49] <= 5
 R225: 0.6 control[11] + control[30] - mode[49] >= -1.183
 R226: 0.6 control[11] + control[30] + mode[49] <= 0.817
 R227: mode[11] + mode[30] + mode[49] = 1
 R228: state[12] - state[13] + 0.1538782696540474 control[12] = 0
 R229: state[32] + 0.0469328722444845 state[52] - state[33]
   + 0.1538782696540474 control[31] = 0
 R230: state[52] - 0.2118055555555556 state[72] - state[53]
   - 0.6944444444444444 control[31] = 0
 R231: state[72] - state[73] + control[50] = 0
 R232: 0.1538782696540474 control[12] >= -0.0469328722444845
 R233: 0.1538782696540474 control[12] <= 0.2530671277555155
 R234: 0.00953125 state[72] + 0.1851282696540474 control[31]
   + 0.045 control[50] >= -0.3
 R235: 0.00953125 state[72] + 0.1851282696540474 control[31]
   + 0.045 control[50] <= 0.3
 R236: - 0.6 control[12] + control[31] <= 0.183
 R237: 0.6 control[12] + control[31] >= -0.183
 R238: control[50] - mode[12] >= -1
 R239: control[50] + mode[12] <= 1
 R240: control[50] - 5 mode[31] >= -5
 R241: 0.6 control[12] - control[31] - mode[31] >= -1.183
 R242: 0.6 control[12] - control[31] + mode[31] <= 0.817
 R243: control[50] + 5 mode[50] <= 5
 R244: 0.6 control[12] + control[31] - mode[50] >= -1.183
 R245: 0.6 control[12] + control[31] + mode[50] <= 0.817
 R246: mode[12] + mode[31] + mode[50] = 1
 R247: state[13] - state[14] + 0.1538782696540474 control[13] = 0
 R248: state[33] + 0.0469328722444845 state[53] - state[34]
   + 0.1538782696540474 control[32] = 0
 R249: state[53] - 0.2118055555555556 state[73] - state[54]
   - 0.6944444444444444 control[32] = 0
 R250: state[73] - state[74] + control[51] = 0
 R251: 0.1538782696540474 control[13] >= -0.0469328722444845
 R252: 0.1538782696540474 control[13] <= 0.2530671277555155
 R253: 0.00953125 state[73] + 0.1851282696540474 control[32]
   + 0.045 control[51] >= -0.3
 R254: 0.00953125 state[73] + 0.1851282696540474 control[32]
   + 0.045 control[51] <= 0.3
 R255: - 0.6 control[13] + control[32] <= 0.183
 R256: 0.6 control[13] + control[32] >= -0.183
 R257: control[51] - mode[13] >= -1
 R258: control[51] + mode[13] <= 1
 R259: control[51] - 5 mode[32] >= -5
 R260: 0.6 control[13] - control[32] - mode[32] >= -1.183
 R261: 0.6 control[13] - control[32] + mode[32] <= 0.817
 R262: control[51] + 5 mode[51] <= 5
 R263: 0.6 control[13] + control[32] - mode[51] >= -1.183
 R264: 0.6 control[13] + control[32] + mode[51] <= 0.817
 R265: mode[13] + mode[32] + mode[51] = 1
 R266: state[14] - state[15] + 0.1538782696540474 control[14] = 0
 R267: state[34] + 0.0469328722444845 state[54] - state[35]
   + 0.1538782696540474 control[33] = 0
 R268: state[54] - 0.2118055555555556 state[74] - state[55]
   - 0.6944444444444444 control[33] = 0
 R269: state[74] - state[75] + control[52] = 0
 R270: 0.1538782696540474 control[14] >= -0.0469328722444845
 R271: 0.1538782696540474 control[14] <= 0.2530671277555155
 R272: 0.00953125 state[74] + 0.1851282696540474 control[33]
   + 0.045 control[52] >= -0.3
 R273: 0.00953125 state[74] + 0.1851282696540474 control[33]
   + 0.045 control[52] <= 0.3
 R274: - 0.6 control[14] + control[33] <= 0.183
 R275: 0.6 control[14] + control[33] >= -0.183
 R276: control[52] - mode[14] >= -1
 R277: control[52] + mode[14] <= 1
 R278: control[52] - 5 mode[33] >= -5
 R279: 0.6 control[14] - control[33] - mode[33] >= -1.183
 R280: 0.6 control[14] - control[33] + mode[33] <= 0.817
 R281: control[52] + 5 mode[52] <= 5
 R282: 0.6 control[14] + control[33] - mode[52] >= -1.183
 R283: 0.6 control[14] + control[33] + mode[52] <= 0.817
 R284: mode[14] + mode[33] + mode[52] = 1
 R285: state[15] - state[16] + 0.1538782696540474 control[15] = 0
 R286: state[35] + 0.0469328722444845 state[55] - state[36]
   + 0.1538782696540474 control[34] = 0
 R287: state[55] - 0.2118055555555556 state[75] - state[56]
   - 0.6944444444444444 control[34] = 0
 R288: state[75] - state[76] + control[53] = 0
 R289: 0.1538782696540474 control[15] >= -0.0469328722444845
 R290: 0.1538782696540474 control[15] <= 0.2530671277555155
 R291: 0.00953125 state[75] + 0.1851282696540474 control[34]
   + 0.045 control[53] >= -0.3
 R292: 0.00953125 state[75] + 0.1851282696540474 control[34]
   + 0.045 control[53] <= 0.3
 R293: - 0.6 control[15] + control[34] <= 0.183
 R294: 0.6 control[15] + control[34] >= -0.183
 R295: control[53] - mode[15] >= -1
 R296: control[53] + mode[15] <= 1
 R297: control[53] - 5 mode[34] >= -5
 R298: 0.6 control[15] - control[34] - mode[34] >= -1.183
 R299: 0.6 control[15] - control[34] + mode[34] <= 0.817
 R300: control[53] + 5 mode[53] <= 5
 R301: 0.6 control[15] + control[34] - mode[53] >= -1.183
 R302: 0.6 control[15] + control[34] + mode[53] <= 0.817
 R303: mode[15] + mode[34] + mode[53] = 1
 R304: state[16] - state[17] + 0.1538782696540474 control[16] = 0
 R305: state[36] + 0.0469328722444845 state[56] - state[37]
   + 0.1538782696540474 control[35] = 0
 R306: state[56] - 0.2118055555555556 state[76] - state[57]
   - 0.6944444444444444 control[35] = 0
 R307: state[76] - state[77] + control[54] = 0
 R308: 0.1538782696540474 control[16] >= -0.0469328722444845
 R309: 0.1538782696540474 control[16] <= 0.2530671277555155
 R310: 0.00953125 state[76] + 0.1851282696540474 control[35]
   + 0.045 control[54] >= -0.3
 R311: 0.00953125 state[76] + 0.1851282696540474 control[35]
   + 0.045 control[54] <= 0.3
 R312: - 0.6 control[16] + control[35] <= 0.183
 R313: 0.6 control[16] + control[35] >= -0.183
 R314: control[54] - mode[16] >= -1
 R315: control[54] + mode[16] <= 1
 R316: control[54] - 5 mode[35] >= -5
 R317: 0.6 control[16] - control[35] - mode[35] >= -1.183
 R318: 0.6 control[16] - control[35] + mode[35] <= 0.817
 R319: control[54] + 5 mode[54] <= 5
 R320: 0.6 control[16] + control[35] - mode[54] >= -1.183
 R321: 0.6 control[16] + control[35] + mode[54] <= 0.817
 R322: mode[16] + mode[35] + mode[54] = 1
 R323: state[17] - state[18] + 0.1538782696540474 control[17] = 0
 R324: state[37] + 0.0469328722444845 state[57] - state[38]
   + 0.1538782696540474 control[36] = 0
 R325: state[57] - 0.2118055555555556 state[77] - state[58]
   - 0.6944444444444444 control[36] = 0
 R326: state[77] - state[78] + control[55] = 0
 R327: 0.1538782696540474 control[17] >= -0.0469328722444845
 R328: 0.1538782696540474 control[17] <= 0.2530671277555155
 R329: 0.00953125 state[77] + 0.1851282696540474 control[36]
   + 0.045 control[55] >= -0.3
 R330: 0.00953125 state[77] + 0.1851282696540474 control[36]
   + 0.045 control[55] <= 0.3
 R331: - 0.6 control[17] + control[36] <= 0.183
 R332: 0.6 control[17] + control[36] >= -0.183
 R333: control[55] - mode[17] >= -1
 R334: control[55] + mode[17] <= 1
 R335: control[55] - 5 mode[36] >= -5
 R336: 0.6 control[17] - control[36] - mode[36] >= -1.183
 R337: 0.6 control[17] - control[36] + mode[36] <= 0.817
 R338: control[55] + 5 mode[55] <= 5
 R339: 0.6 control[17] + control[36] - mode[55] >= -1.183
 R340: 0.6 control[17] + control[36] + mode[55] <= 0.817
 R341: mode[17] + mode[36] + mode[55] = 1
 R342: state[18] - state[19] + 0.1538782696540474 control[18] = 0
 R343: state[38] + 0.0469328722444845 state[58] - state[39]
   + 0.1538782696540474 control[37] = 0
 R344: state[58] - 0.2118055555555556 state[78] - state[59]
   - 0.6944444444444444 control[37] = 0
 R345: state[78] - state[79] + control[56] = 0
 R346: 0.1538782696540474 control[18] >= -0.0469328722444845
 R347: 0.1538782696540474 control[18] <= 0.2530671277555155
 R348: 0.00953125 state[78] + 0.1851282696540474 control[37]
   + 0.045 control[56] >= -0.3
 R349: 0.00953125 state[78] + 0.1851282696540474 control[37]
   + 0.045 control[56] <= 0.3
 R350: - 0.6 control[18] + control[37] <= 0.183
 R351: 0.6 control[18] + control[37] >= -0.183
 R352: control[56] - mode[18] >= -1
 R353: control[56] + mode[18] <= 1
 R354: control[56] - 5 mode[37] >= -5
 R355: 0.6 control[18] - control[37] - mode[37] >= -1.183
 R356: 0.6 control[18] - control[37] + mode[37] <= 0.817
 R357: control[56] + 5 mode[56] <= 5
 R358: 0.6 control[18] + control[37] - mode[56] >= -1.183
 R359: 0.6 control[18] + control[37] + mode[56] <= 0.817
 R360: mode[18] + mode[37] + mode[56] = 1
Bounds
 state[0] = 0.00394623403482919
 state[20] = 0.0242826197296378
 state[40] = 0.00843130052089691
 state[60] = 5.4069316138275175e-04
 state[1] free
 state[21] free
 state[41] free
 -0.7 <= state[61] <= 0.7
 state[2] free
 state[22] free
 state[42] free
 -0.7 <= state[62] <= 0.7
 state[3] free
 state[23] free
 state[43] free
 -0.7 <= state[63] <= 0.7
 state[4] free
 state[24] free
 state[44] free
 -0.7 <= state[64] <= 0.7
 state[5] free
 state[25] free
 state[45] free
 -0.7 <= state[65] <= 0.7
 state[6] free
 state[26] free
 state[46] free
 -0.7 <= state[66] <= 0.7
 state[7] free
 state[27] free
 state[47] free
 -0.7 <= state[67] <= 0.7
 state[8] free
 state[28] free
 state[48] free
 -0.7 <= state[68] <= 0.7
 state[9] free
 state[29] free
 state[49] free
 -0.7 <= state[69] <= 0.7
 state[10] free
 state[30] free
 state[50] free
 -0.7 <= state[70] <= 0.7
 state[11] free
 state[31] free
 state[51] free
 -0.7 <= state[71] <= 0.7
 state[12] free
 state[32] free
 state[52] free
 -0.7 <= state[72] <= 0.7
 state[13] free
 state[33] free
 state[53] free
 -0.7 <= state[73] <= 0.7
 state[14] free
 state[34] free
 state[54] free
 -0.7 <= state[74] <= 0.7
 state[15] free
 state[35] free
 state[55] free
 -0.7 <= state[75] <= 0.7
 state[16] free
 state[36] free
 state[56] free
 -0.7 <= state[76] <= 0.7
 state[17] free
 state[37] free
 state[57] free
 -0.7 <= state[77] <= 0.7
 state[18] free
 state[38] free
 state[58] free
 -0.7 <= state[78] <= 0.7
 state[19] free
 state[39] free
 state[59] free
 -0.7 <= state[79] <= 0.7
 control[0] >= -0.305
 control[19] free
 control[38] free
 control[1] >= -0.305
 control[20] free
 control[39] free
 control[2] >= -0.305
 control[21] free
 control[40] free
 control[3] >= -0.305
 control[22] free
 control[41] free
 control[4] >= -0.305
 control[23] free
 control[42] free
 control[5] >= -0.305
 control[24] free
 control[43] free
 control[6] >= -0.305
 control[25] free
 control[44] free
 control[7] >= -0.305
 control[26] free
 control[45] free
 control[8] >= -0.305
 control[27] free
 control[46] free
 control[9] >= -0.305
 control[28] free
 control[47] free
 control[10] >= -0.305
 control[29] free
 control[48] free
 control[11] >= -0.305
 control[30] free
 control[49] free
 control[12] >= -0.305
 control[31] free
 control[50] free
 control[13] >= -0.305
 control[32] free
 control[51] free
 control[14] >= -0.305
 control[33] free
 control[52] free
 control[15] >= -0.305
 control[34] free
 control[53] free
 control[16] >= -0.305
 control[35] free
 control[54] free
 control[17] >= -0.305
 control[36] free
 control[55] free
 control[18] >= -0.305
 control[37] free
 control[56] free
Binaries
 mode[0] mode[19] mode[38] mode[1] mode[20] mode[39] mode[2] mode[21]
 mode[40] mode[3] mode[22] mode[41] mode[4] mode[23] mode[42] mode[5]
 mode[24] mode[43] mode[6] mode[25] mode[44] mode[7] mode[26] mode[45]
 mode[8] mode[27] mode[46] mode[9] mode[28] mode[47] mode[10] mode[29]
 mode[48] mode[11] mode[30] mode[49] mode[12] mode[31] mode[50] mode[13]
 mode[32] mode[51] mode[14] mode[33] mode[52] mode[15] mode[34] mode[53]
 mode[16] mode[35] mode[54] mode[17] mode[36] mode[55] mode[18] mode[37]
 mode[56]
End
