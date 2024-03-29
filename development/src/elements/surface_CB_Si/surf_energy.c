#include "Tersoff_inc_surf.h"

#include <math.h>

static double z[437];

/* calculate the surface strain energy density */
double get_energy_surf(const double* params, const double *Xsi, const double *Xa, const double *Ya, const double *Za, const double* Cmat) {

/* common definitions */
#include "Tersoff_common_defines_surf.h"

	z[1] = c*c;
	z[2] = 1./(d*d);
	z[3] = d*d;
	z[4] = 1./n;
	z[5] = -X1;
	z[6] = -X2;
	z[7] = -X3;
	z[8] = -X4;
	z[9] = -X5;
	z[10] = -X6;
	z[11] = -X7;
	z[12] = -X8;
	z[13] = -X9;
	z[14] = -Xs1;
	z[15] = -Y1;
	z[16] = -Y2;
	z[17] = -Y3;
	z[18] = -Y4;
	z[19] = -Y5;
	z[20] = -Y6;
	z[21] = -Y7;
	z[22] = -Y8;
	z[23] = -Y9;
	z[24] = -Ys1;
	z[25] = -Z1;
	z[26] = -Z2;
	z[27] = -Z3;
	z[28] = -Z4;
	z[29] = -Z5;
	z[30] = -Z6;
	z[31] = -Z7;
	z[32] = -Z8;
	z[33] = -Z9;
	z[34] = -Zs1;
	z[35] = X2 + z[10];
	z[10] = X1 + Xs1 + z[10];
	z[36] = X1 + z[11];
	z[37] = X8 + z[11];
	z[38] = X9 + z[11];
	z[39] = X1 + z[12];
	z[40] = X7 + z[12];
	z[41] = X9 + z[12];
	z[42] = X1 + z[13];
	z[43] = X7 + z[13];
	z[44] = X8 + z[13];
	z[11] = X6 + z[11] + z[14];
	z[12] = X6 + z[12] + z[14];
	z[13] = X6 + z[13] + z[14];
	z[45] = Y3 + z[15];
	z[46] = Y4 + z[15];
	z[47] = Y5 + z[15];
	z[48] = Y7 + z[15];
	z[49] = Y8 + z[15];
	z[50] = Y9 + z[15];
	z[51] = Y6 + z[16];
	z[16] = Y1 + Ys1 + z[16];
	z[52] = Y1 + z[17];
	z[53] = Y4 + z[17];
	z[54] = Y5 + z[17];
	z[55] = Y1 + z[18];
	z[56] = Y3 + z[18];
	z[57] = Y5 + z[18];
	z[58] = Y1 + z[19];
	z[59] = Y3 + z[19];
	z[60] = Y4 + z[19];
	z[2] = z[1]*z[2];
	z[61] = Y2 + z[20];
	z[20] = Y1 + Ys1 + z[20];
	z[62] = Y1 + z[21];
	z[63] = Y8 + z[21];
	z[64] = Y9 + z[21];
	z[65] = Y1 + z[22];
	z[66] = Y7 + z[22];
	z[67] = Y9 + z[22];
	z[68] = Y1 + z[23];
	z[69] = Y7 + z[23];
	z[70] = Y8 + z[23];
	z[71] = Y2 + z[15] + z[24];
	z[15] = Y6 + z[15] + z[24];
	z[17] = Y2 + z[17] + z[24];
	z[18] = Y2 + z[18] + z[24];
	z[19] = Y2 + z[19] + z[24];
	z[21] = Y6 + z[21] + z[24];
	z[22] = Y6 + z[22] + z[24];
	z[23] = Y6 + z[23] + z[24];
	z[24] = Z1 + z[27];
	z[72] = Z1 + z[28];
	z[73] = Z1 + z[29];
	z[74] = z[25] + Z3;
	z[75] = z[28] + Z3;
	z[76] = z[29] + Z3;
	z[77] = Z2 + z[30];
	z[78] = Z1 + z[31];
	z[79] = Z1 + z[32];
	z[80] = Z1 + z[33];
	z[81] = Z2 + z[25] + z[34];
	z[82] = Z2 + z[27] + z[34];
	z[83] = Z2 + z[28] + z[34];
	z[84] = Z2 + z[29] + z[34];
	z[4] = -0.5*z[4];
	z[85] = z[25] + Z4;
	z[86] = z[27] + Z4;
	z[29] = z[29] + Z4;
	z[87] = X3 + z[5];
	z[88] = X4 + z[5];
	z[89] = X5 + z[5];
	z[90] = X7 + z[5];
	z[91] = X8 + z[5];
	z[92] = X9 + z[5];
	z[93] = X2 + z[14] + z[5];
	z[5] = X6 + z[14] + z[5];
	z[94] = z[25] + Z5;
	z[27] = z[27] + Z5;
	z[28] = z[28] + Z5;
	z[95] = X6 + z[6];
	z[6] = X1 + Xs1 + z[6];
	z[96] = z[26] + Z6;
	z[97] = z[25] + z[34] + Z6;
	z[98] = z[31] + z[34] + Z6;
	z[99] = z[32] + z[34] + Z6;
	z[34] = z[33] + z[34] + Z6;
	z[100] = X1 + z[7];
	z[101] = X4 + z[7];
	z[102] = X5 + z[7];
	z[7] = X2 + z[14] + z[7];
	z[103] = z[25] + Z7;
	z[104] = z[32] + Z7;
	z[105] = z[33] + Z7;
	z[106] = X1 + z[8];
	z[107] = X3 + z[8];
	z[108] = X5 + z[8];
	z[8] = X2 + z[14] + z[8];
	z[109] = z[25] + Z8;
	z[110] = z[31] + Z8;
	z[33] = z[33] + Z8;
	z[111] = X1 + z[9];
	z[112] = X3 + z[9];
	z[113] = X4 + z[9];
	z[9] = X2 + z[14] + z[9];
	z[14] = z[25] + Z9;
	z[25] = z[31] + Z9;
	z[31] = z[32] + Z9;
	z[26] = Z1 + z[26] + Zs1;
	z[30] = Z1 + z[30] + Zs1;
	z[32] = C31*z[85];
	z[114] = C32*z[85];
	z[115] = C33*z[85];
	z[116] = C31*z[86];
	z[117] = C32*z[86];
	z[118] = C33*z[86];
	z[119] = C31*z[29];
	z[120] = C32*z[29];
	z[121] = C33*z[29];
	z[122] = C11*z[87];
	z[123] = C12*z[87];
	z[124] = C13*z[87];
	z[125] = C11*z[88];
	z[126] = C12*z[88];
	z[127] = C13*z[88];
	z[128] = C11*z[89];
	z[129] = C12*z[89];
	z[130] = C13*z[89];
	z[131] = C11*z[90];
	z[132] = C12*z[90];
	z[133] = C13*z[90];
	z[134] = C11*z[91];
	z[135] = C12*z[91];
	z[136] = C13*z[91];
	z[137] = C11*z[92];
	z[138] = C12*z[92];
	z[139] = C13*z[92];
	z[140] = C11*z[93];
	z[141] = C12*z[93];
	z[142] = C13*z[93];
	z[143] = C11*z[5];
	z[144] = C12*z[5];
	z[145] = C13*z[5];
	z[146] = C31*z[94];
	z[147] = C32*z[94];
	z[148] = C33*z[94];
	z[149] = C31*z[27];
	z[150] = C32*z[27];
	z[151] = C33*z[27];
	z[152] = C31*z[28];
	z[153] = C32*z[28];
	z[154] = C33*z[28];
	z[155] = C11*z[95];
	z[156] = C12*z[95];
	z[157] = C13*z[95];
	z[158] = C11*z[6];
	z[159] = C12*z[6];
	z[160] = C13*z[6];
	z[161] = C31*z[96];
	z[162] = C32*z[96];
	z[163] = C33*z[96];
	z[164] = C31*z[97];
	z[165] = C32*z[97];
	z[166] = C33*z[97];
	z[167] = C31*z[98];
	z[168] = C32*z[98];
	z[169] = C33*z[98];
	z[170] = C31*z[99];
	z[171] = C32*z[99];
	z[172] = C33*z[99];
	z[173] = C31*z[34];
	z[174] = C32*z[34];
	z[175] = C33*z[34];
	z[176] = C11*z[100];
	z[177] = C12*z[100];
	z[178] = C13*z[100];
	z[179] = C11*z[101];
	z[180] = C12*z[101];
	z[181] = C13*z[101];
	z[182] = C11*z[102];
	z[183] = C12*z[102];
	z[184] = C13*z[102];
	z[185] = C11*z[7];
	z[186] = C12*z[7];
	z[187] = C13*z[7];
	z[188] = C31*z[103];
	z[189] = C32*z[103];
	z[190] = C33*z[103];
	z[191] = C31*z[104];
	z[192] = C32*z[104];
	z[193] = C33*z[104];
	z[194] = C31*z[105];
	z[195] = C32*z[105];
	z[196] = C33*z[105];
	z[197] = C11*z[106];
	z[198] = C12*z[106];
	z[199] = C13*z[106];
	z[200] = C11*z[107];
	z[201] = C12*z[107];
	z[202] = C13*z[107];
	z[203] = C11*z[108];
	z[204] = C12*z[108];
	z[205] = C13*z[108];
	z[206] = C11*z[8];
	z[207] = C12*z[8];
	z[208] = C13*z[8];
	z[209] = C31*z[109];
	z[210] = C32*z[109];
	z[211] = C33*z[109];
	z[212] = C31*z[110];
	z[213] = C32*z[110];
	z[214] = C33*z[110];
	z[215] = C31*z[33];
	z[216] = C32*z[33];
	z[217] = C33*z[33];
	z[218] = C11*z[111];
	z[219] = C12*z[111];
	z[220] = C13*z[111];
	z[221] = C11*z[112];
	z[222] = C12*z[112];
	z[223] = C13*z[112];
	z[224] = C11*z[113];
	z[225] = C12*z[113];
	z[226] = C13*z[113];
	z[227] = C11*z[9];
	z[228] = C12*z[9];
	z[229] = C13*z[9];
	z[230] = C31*z[14];
	z[231] = C32*z[14];
	z[232] = C33*z[14];
	z[233] = C31*z[25];
	z[234] = C32*z[25];
	z[235] = C33*z[25];
	z[236] = C31*z[31];
	z[237] = C32*z[31];
	z[238] = C33*z[31];
	z[239] = C31*z[26];
	z[240] = C32*z[26];
	z[241] = C33*z[26];
	z[242] = C31*z[30];
	z[243] = C32*z[30];
	z[244] = C33*z[30];
	z[245] = C11*z[35];
	z[246] = C12*z[35];
	z[247] = C13*z[35];
	z[248] = C11*z[10];
	z[249] = C12*z[10];
	z[250] = C13*z[10];
	z[251] = C11*z[36];
	z[252] = C12*z[36];
	z[253] = C13*z[36];
	z[254] = C11*z[37];
	z[255] = C12*z[37];
	z[256] = C13*z[37];
	z[257] = C11*z[38];
	z[258] = C12*z[38];
	z[259] = C13*z[38];
	z[260] = C11*z[39];
	z[261] = C12*z[39];
	z[262] = C13*z[39];
	z[263] = C11*z[40];
	z[264] = C12*z[40];
	z[265] = C13*z[40];
	z[266] = C11*z[41];
	z[267] = C12*z[41];
	z[268] = C13*z[41];
	z[269] = C11*z[42];
	z[270] = C12*z[42];
	z[271] = C13*z[42];
	z[272] = C11*z[43];
	z[273] = C12*z[43];
	z[274] = C13*z[43];
	z[275] = C11*z[44];
	z[276] = C12*z[44];
	z[277] = C13*z[44];
	z[278] = C11*z[11];
	z[279] = C12*z[11];
	z[280] = C13*z[11];
	z[281] = C11*z[12];
	z[282] = C12*z[12];
	z[283] = C13*z[12];
	z[284] = C11*z[13];
	z[285] = C12*z[13];
	z[286] = C13*z[13];
	z[287] = C21*z[45];
	z[288] = C22*z[45];
	z[289] = C23*z[45];
	z[290] = C21*z[46];
	z[291] = C22*z[46];
	z[292] = C23*z[46];
	z[293] = C21*z[47];
	z[294] = C22*z[47];
	z[295] = C23*z[47];
	z[296] = C21*z[48];
	z[297] = C22*z[48];
	z[298] = C23*z[48];
	z[299] = C21*z[49];
	z[300] = C22*z[49];
	z[301] = C23*z[49];
	z[302] = C21*z[50];
	z[303] = C22*z[50];
	z[304] = C23*z[50];
	z[305] = C21*z[51];
	z[306] = C22*z[51];
	z[307] = C23*z[51];
	z[308] = C21*z[16];
	z[309] = C22*z[16];
	z[310] = C23*z[16];
	z[311] = C21*z[52];
	z[312] = C22*z[52];
	z[313] = C23*z[52];
	z[314] = C21*z[53];
	z[315] = C22*z[53];
	z[316] = C23*z[53];
	z[317] = C21*z[54];
	z[318] = C22*z[54];
	z[319] = C23*z[54];
	z[320] = C21*z[55];
	z[321] = C22*z[55];
	z[322] = C23*z[55];
	z[323] = C21*z[56];
	z[324] = C22*z[56];
	z[325] = C23*z[56];
	z[326] = C21*z[57];
	z[327] = C22*z[57];
	z[328] = C23*z[57];
	z[329] = C21*z[58];
	z[330] = C22*z[58];
	z[331] = C23*z[58];
	z[332] = C21*z[59];
	z[333] = C22*z[59];
	z[334] = C23*z[59];
	z[335] = C21*z[60];
	z[336] = C22*z[60];
	z[337] = C23*z[60];
	z[338] = C21*z[61];
	z[339] = C22*z[61];
	z[340] = C23*z[61];
	z[341] = C21*z[20];
	z[342] = C22*z[20];
	z[343] = C23*z[20];
	z[344] = C21*z[62];
	z[345] = C22*z[62];
	z[346] = C23*z[62];
	z[347] = C21*z[63];
	z[348] = C22*z[63];
	z[349] = C23*z[63];
	z[350] = C21*z[64];
	z[351] = C22*z[64];
	z[352] = C23*z[64];
	z[353] = C21*z[65];
	z[354] = C22*z[65];
	z[355] = C23*z[65];
	z[356] = C21*z[66];
	z[357] = C22*z[66];
	z[358] = C23*z[66];
	z[359] = C21*z[67];
	z[360] = C22*z[67];
	z[361] = C23*z[67];
	z[362] = C21*z[68];
	z[363] = C22*z[68];
	z[364] = C23*z[68];
	z[365] = C21*z[69];
	z[366] = C22*z[69];
	z[367] = C23*z[69];
	z[368] = C21*z[70];
	z[369] = C22*z[70];
	z[370] = C23*z[70];
	z[371] = C21*z[71];
	z[372] = C22*z[71];
	z[373] = C23*z[71];
	z[374] = C21*z[15];
	z[375] = C22*z[15];
	z[376] = C23*z[15];
	z[377] = C21*z[17];
	z[378] = C22*z[17];
	z[379] = C23*z[17];
	z[380] = C21*z[18];
	z[381] = C22*z[18];
	z[382] = C23*z[18];
	z[383] = C21*z[19];
	z[384] = C22*z[19];
	z[385] = C23*z[19];
	z[386] = C21*z[21];
	z[387] = C22*z[21];
	z[388] = C23*z[21];
	z[389] = C21*z[22];
	z[390] = C22*z[22];
	z[391] = C23*z[22];
	z[392] = C21*z[23];
	z[393] = C22*z[23];
	z[394] = C23*z[23];
	z[395] = C31*z[24];
	z[396] = C32*z[24];
	z[397] = C33*z[24];
	z[398] = C31*z[72];
	z[399] = C32*z[72];
	z[400] = C33*z[72];
	z[401] = C31*z[73];
	z[402] = C32*z[73];
	z[403] = C33*z[73];
	z[404] = C31*z[74];
	z[405] = C32*z[74];
	z[406] = C33*z[74];
	z[407] = C31*z[75];
	z[408] = C32*z[75];
	z[409] = C33*z[75];
	z[410] = C31*z[76];
	z[411] = C32*z[76];
	z[412] = C33*z[76];
	z[413] = C31*z[77];
	z[414] = C32*z[77];
	z[415] = C33*z[77];
	z[416] = C31*z[78];
	z[417] = C32*z[78];
	z[418] = C33*z[78];
	z[419] = C31*z[79];
	z[420] = C32*z[79];
	z[421] = C33*z[79];
	z[422] = C31*z[80];
	z[423] = C32*z[80];
	z[424] = C33*z[80];
	z[425] = C31*z[81];
	z[426] = C32*z[81];
	z[427] = C33*z[81];
	z[428] = C31*z[82];
	z[429] = C32*z[82];
	z[430] = C33*z[82];
	z[431] = C31*z[83];
	z[432] = C32*z[83];
	z[433] = C33*z[83];
	z[434] = C31*z[84];
	z[435] = C32*z[84];
	z[436] = C33*z[84];
	z[32] = z[125] + z[290] + z[32];
	z[114] = z[114] + z[126] + z[291];
	z[115] = z[115] + z[127] + z[292];
	z[125] = z[128] + z[146] + z[293];
	z[126] = z[129] + z[147] + z[294];
	z[127] = z[130] + z[148] + z[295];
	z[128] = z[131] + z[188] + z[296];
	z[129] = z[132] + z[189] + z[297];
	z[130] = z[133] + z[190] + z[298];
	z[131] = z[134] + z[209] + z[299];
	z[132] = z[135] + z[210] + z[300];
	z[133] = z[136] + z[211] + z[301];
	z[134] = z[137] + z[230] + z[302];
	z[135] = z[138] + z[231] + z[303];
	z[136] = z[139] + z[232] + z[304];
	z[137] = z[155] + z[161] + z[305];
	z[138] = z[156] + z[162] + z[306];
	z[139] = z[157] + z[163] + z[307];
	z[146] = z[158] + z[239] + z[308];
	z[147] = z[159] + z[240] + z[309];
	z[148] = z[160] + z[241] + z[310];
	z[116] = z[116] + z[179] + z[314];
	z[117] = z[117] + z[180] + z[315];
	z[118] = z[118] + z[181] + z[316];
	z[149] = z[149] + z[182] + z[317];
	z[150] = z[150] + z[183] + z[318];
	z[151] = z[151] + z[184] + z[319];
	z[152] = z[152] + z[203] + z[326];
	z[153] = z[153] + z[204] + z[327];
	z[154] = z[154] + z[205] + z[328];
	z[119] = z[119] + z[224] + z[335];
	z[120] = z[120] + z[225] + z[336];
	z[121] = z[121] + z[226] + z[337];
	z[155] = z[242] + z[248] + z[341];
	z[156] = z[243] + z[249] + z[342];
	z[157] = z[244] + z[250] + z[343];
	z[158] = z[212] + z[254] + z[347];
	z[159] = z[213] + z[255] + z[348];
	z[160] = z[214] + z[256] + z[349];
	z[161] = z[233] + z[257] + z[350];
	z[162] = z[234] + z[258] + z[351];
	z[163] = z[235] + z[259] + z[352];
	z[179] = z[191] + z[263] + z[356];
	z[180] = z[192] + z[264] + z[357];
	z[181] = z[193] + z[265] + z[358];
	z[182] = z[236] + z[266] + z[359];
	z[183] = z[237] + z[267] + z[360];
	z[184] = z[238] + z[268] + z[361];
	z[188] = z[194] + z[272] + z[365];
	z[189] = z[195] + z[273] + z[366];
	z[190] = z[196] + z[274] + z[367];
	z[191] = z[215] + z[275] + z[368];
	z[192] = z[216] + z[276] + z[369];
	z[193] = z[217] + z[277] + z[370];
	z[143] = z[143] + z[164] + z[374];
	z[144] = z[144] + z[165] + z[375];
	z[145] = z[145] + z[166] + z[376];
	z[164] = z[167] + z[278] + z[386];
	z[165] = z[168] + z[279] + z[387];
	z[166] = z[169] + z[280] + z[388];
	z[167] = z[170] + z[281] + z[389];
	z[168] = z[171] + z[282] + z[390];
	z[169] = z[172] + z[283] + z[391];
	z[170] = z[173] + z[284] + z[392];
	z[171] = z[174] + z[285] + z[393];
	z[172] = z[175] + z[286] + z[394];
	z[173] = z[176] + z[311] + z[395];
	z[174] = z[177] + z[312] + z[396];
	z[175] = z[178] + z[313] + z[397];
	z[176] = z[197] + z[320] + z[398];
	z[177] = z[198] + z[321] + z[399];
	z[178] = z[199] + z[322] + z[400];
	z[194] = z[218] + z[329] + z[401];
	z[195] = z[219] + z[330] + z[402];
	z[196] = z[220] + z[331] + z[403];
	z[122] = z[122] + z[287] + z[404];
	z[123] = z[123] + z[288] + z[405];
	z[124] = z[124] + z[289] + z[406];
	z[197] = z[200] + z[323] + z[407];
	z[198] = z[201] + z[324] + z[408];
	z[199] = z[202] + z[325] + z[409];
	z[200] = z[221] + z[332] + z[410];
	z[201] = z[222] + z[333] + z[411];
	z[202] = z[223] + z[334] + z[412];
	z[203] = z[245] + z[338] + z[413];
	z[204] = z[246] + z[339] + z[414];
	z[205] = z[247] + z[340] + z[415];
	z[209] = z[251] + z[344] + z[416];
	z[210] = z[252] + z[345] + z[417];
	z[211] = z[253] + z[346] + z[418];
	z[212] = z[260] + z[353] + z[419];
	z[213] = z[261] + z[354] + z[420];
	z[214] = z[262] + z[355] + z[421];
	z[215] = z[269] + z[362] + z[422];
	z[216] = z[270] + z[363] + z[423];
	z[217] = z[271] + z[364] + z[424];
	z[140] = z[140] + z[371] + z[425];
	z[141] = z[141] + z[372] + z[426];
	z[142] = z[142] + z[373] + z[427];
	z[185] = z[185] + z[377] + z[428];
	z[186] = z[186] + z[378] + z[429];
	z[187] = z[187] + z[379] + z[430];
	z[206] = z[206] + z[380] + z[431];
	z[207] = z[207] + z[381] + z[432];
	z[208] = z[208] + z[382] + z[433];
	z[218] = z[227] + z[383] + z[434];
	z[219] = z[228] + z[384] + z[435];
	z[220] = z[229] + z[385] + z[436];
	z[32] = -z[32]*z[88];
	z[85] = -z[115]*z[85];
	z[88] = -z[125]*z[89];
	z[89] = -z[127]*z[94];
	z[90] = -z[128]*z[90];
	z[94] = -z[103]*z[130];
	z[91] = -z[131]*z[91];
	z[103] = -z[109]*z[133];
	z[92] = -z[134]*z[92];
	z[14] = -z[136]*z[14];
	z[95] = -z[137]*z[95];
	z[96] = -z[139]*z[96];
	z[6] = z[146]*z[6];
	z[26] = z[148]*z[26];
	z[101] = -z[101]*z[116];
	z[86] = -z[118]*z[86];
	z[102] = -z[102]*z[149];
	z[27] = -z[151]*z[27];
	z[108] = -z[108]*z[152];
	z[28] = -z[154]*z[28];
	z[109] = -z[113]*z[119];
	z[46] = -z[114]*z[46];
	z[29] = -z[121]*z[29];
	z[10] = z[10]*z[155];
	z[30] = z[157]*z[30];
	z[37] = -z[158]*z[37];
	z[110] = -z[110]*z[160];
	z[38] = -z[161]*z[38];
	z[47] = -z[126]*z[47];
	z[25] = -z[163]*z[25];
	z[40] = -z[179]*z[40];
	z[104] = -z[104]*z[181];
	z[41] = -z[182]*z[41];
	z[31] = -z[184]*z[31];
	z[43] = -z[188]*z[43];
	z[105] = -z[105]*z[190];
	z[48] = -z[129]*z[48];
	z[44] = -z[191]*z[44];
	z[33] = -z[193]*z[33];
	z[5] = z[143]*z[5];
	z[97] = z[145]*z[97];
	z[11] = z[11]*z[164];
	z[98] = z[166]*z[98];
	z[12] = z[12]*z[167];
	z[49] = -z[132]*z[49];
	z[99] = z[169]*z[99];
	z[13] = z[13]*z[170];
	z[34] = z[172]*z[34];
	z[100] = -z[100]*z[173];
	z[106] = -z[106]*z[176];
	z[50] = -z[135]*z[50];
	z[111] = -z[111]*z[194];
	z[87] = -z[122]*z[87];
	z[45] = -z[123]*z[45];
	z[107] = -z[107]*z[197];
	z[51] = -z[138]*z[51];
	z[112] = -z[112]*z[200];
	z[35] = -z[203]*z[35];
	z[36] = -z[209]*z[36];
	z[39] = -z[212]*z[39];
	z[16] = z[147]*z[16];
	z[42] = -z[215]*z[42];
	z[93] = z[140]*z[93];
	z[7] = z[185]*z[7];
	z[52] = -z[174]*z[52];
	z[8] = z[206]*z[8];
	z[9] = z[218]*z[9];
	z[53] = -z[117]*z[53];
	z[54] = -z[150]*z[54];
	z[55] = -z[177]*z[55];
	z[56] = -z[198]*z[56];
	z[57] = -z[153]*z[57];
	z[58] = -z[195]*z[58];
	z[59] = -z[201]*z[59];
	z[60] = -z[120]*z[60];
	z[61] = -z[204]*z[61];
	z[20] = z[156]*z[20];
	z[62] = -z[210]*z[62];
	z[63] = -z[159]*z[63];
	z[64] = -z[162]*z[64];
	z[65] = -z[213]*z[65];
	z[66] = -z[180]*z[66];
	z[67] = -z[183]*z[67];
	z[68] = -z[216]*z[68];
	z[69] = -z[189]*z[69];
	z[70] = -z[192]*z[70];
	z[71] = z[141]*z[71];
	z[15] = z[144]*z[15];
	z[17] = z[17]*z[186];
	z[18] = z[18]*z[207];
	z[19] = z[19]*z[219];
	z[21] = z[165]*z[21];
	z[22] = z[168]*z[22];
	z[23] = z[171]*z[23];
	z[24] = -z[175]*z[24];
	z[72] = -z[178]*z[72];
	z[73] = -z[196]*z[73];
	z[74] = -z[124]*z[74];
	z[75] = -z[199]*z[75];
	z[76] = -z[202]*z[76];
	z[77] = -z[205]*z[77];
	z[78] = -z[211]*z[78];
	z[79] = -z[214]*z[79];
	z[80] = -z[217]*z[80];
	z[81] = z[142]*z[81];
	z[82] = z[187]*z[82];
	z[83] = z[208]*z[83];
	z[84] = z[220]*z[84];
	z[6] = z[16] + z[26] + z[6];
	z[10] = z[10] + z[20] + z[30];
	z[5] = z[15] + z[5] + z[97];
	z[11] = z[11] + z[21] + z[98];
	z[12] = z[12] + z[22] + z[99];
	z[13] = z[13] + z[23] + z[34];
	z[15] = z[71] + z[81] + z[93];
	z[7] = z[17] + z[7] + z[82];
	z[8] = z[18] + z[8] + z[83];
	z[9] = z[19] + z[84] + z[9];
	z[16] = z[10] + z[51] + z[6] + z[95] + z[96];
	z[17] = z[10] + z[35] + z[6] + z[61] + z[77];
	z[18] = z[11] + z[48] + z[5] + z[90] + z[94];
	z[19] = z[11] + z[36] + z[5] + z[62] + z[78];
	z[20] = z[103] + z[12] + z[49] + z[5] + z[91];
	z[21] = z[12] + z[39] + z[5] + z[65] + z[79];
	z[22] = z[11] + z[110] + z[12] + z[37] + z[63];
	z[23] = z[104] + z[11] + z[12] + z[40] + z[66];
	z[14] = z[13] + z[14] + z[5] + z[50] + z[92];
	z[26] = z[13] + z[42] + z[5] + z[68] + z[80];
	z[25] = z[11] + z[13] + z[25] + z[38] + z[64];
	z[30] = z[105] + z[11] + z[13] + z[43] + z[69];
	z[31] = z[12] + z[13] + z[31] + z[41] + z[67];
	z[33] = z[12] + z[13] + z[33] + z[44] + z[70];
	z[24] = z[100] + z[15] + z[24] + z[52] + z[7];
	z[34] = z[15] + z[45] + z[7] + z[74] + z[87];
	z[32] = z[15] + z[32] + z[46] + z[8] + z[85];
	z[35] = z[106] + z[15] + z[55] + z[72] + z[8];
	z[36] = z[101] + z[53] + z[7] + z[8] + z[86];
	z[37] = z[107] + z[56] + z[7] + z[75] + z[8];
	z[38] = z[15] + z[47] + z[88] + z[89] + z[9];
	z[39] = z[111] + z[15] + z[58] + z[73] + z[9];
	z[27] = z[102] + z[27] + z[54] + z[7] + z[9];
	z[40] = z[112] + z[59] + z[7] + z[76] + z[9];
	z[28] = z[108] + z[28] + z[57] + z[8] + z[9];
	z[29] = z[109] + z[29] + z[60] + z[8] + z[9];
	z[41] = 1./sqrt(z[6]);
	z[6] = sqrt(z[6]);
	z[42] = 1./sqrt(z[10]);
	z[10] = sqrt(z[10]);
	z[43] = 1./sqrt(z[5]);
	z[5] = sqrt(z[5]);
	z[44] = 1./sqrt(z[11]);
	z[11] = sqrt(z[11]);
	z[45] = 1./sqrt(z[12]);
	z[12] = sqrt(z[12]);
	z[46] = 1./sqrt(z[13]);
	z[13] = sqrt(z[13]);
	z[47] = 1./sqrt(z[15]);
	z[15] = sqrt(z[15]);
	z[48] = 1./sqrt(z[7]);
	z[7] = sqrt(z[7]);
	z[49] = 1./sqrt(z[8]);
	z[8] = sqrt(z[8]);
	z[50] = 1./sqrt(z[9]);
	z[9] = sqrt(z[9]);
	z[6] = -z[6];
	z[16] = -0.5*z[16]*z[41]*z[42];
	z[17] = -0.5*z[17]*z[41]*z[42];
	z[41] = -lam*z[10];
	z[10] = -mu*z[10];
	z[42] = -lam*z[5];
	z[5] = -mu*z[5];
	z[18] = -0.5*z[18]*z[43]*z[44];
	z[19] = -0.5*z[19]*z[43]*z[44];
	z[51] = -lam*z[11];
	z[11] = -mu*z[11];
	z[20] = -0.5*z[20]*z[43]*z[45];
	z[21] = -0.5*z[21]*z[43]*z[45];
	z[22] = -0.5*z[22]*z[44]*z[45];
	z[23] = -0.5*z[23]*z[44]*z[45];
	z[52] = -lam*z[12];
	z[12] = -mu*z[12];
	z[14] = -0.5*z[14]*z[43]*z[46];
	z[26] = -0.5*z[26]*z[43]*z[46];
	z[25] = -0.5*z[25]*z[44]*z[46];
	z[30] = -0.5*z[30]*z[44]*z[46];
	z[31] = -0.5*z[31]*z[45]*z[46];
	z[33] = -0.5*z[33]*z[45]*z[46];
	z[43] = -lam*z[13];
	z[13] = -mu*z[13];
	z[44] = -lam*z[15];
	z[15] = -mu*z[15];
	z[24] = -0.5*z[24]*z[47]*z[48];
	z[34] = -0.5*z[34]*z[47]*z[48];
	z[45] = -lam*z[7];
	z[7] = -mu*z[7];
	z[32] = -0.5*z[32]*z[47]*z[49];
	z[35] = -0.5*z[35]*z[47]*z[49];
	z[36] = -0.5*z[36]*z[48]*z[49];
	z[37] = -0.5*z[37]*z[48]*z[49];
	z[46] = -lam*z[8];
	z[8] = -mu*z[8];
	z[38] = -0.5*z[38]*z[47]*z[50];
	z[39] = -0.5*z[39]*z[47]*z[50];
	z[27] = -0.5*z[27]*z[48]*z[50];
	z[40] = -0.5*z[40]*z[48]*z[50];
	z[28] = -0.5*z[28]*z[49]*z[50];
	z[29] = -0.5*z[29]*z[49]*z[50];
	z[47] = -lam*z[9];
	z[9] = -mu*z[9];
	z[48] = lam*z[6];
	z[6] = mu*z[6];
	z[41] = exp(z[41]);
	z[10] = exp(z[10]);
	z[42] = exp(z[42]);
	z[5] = exp(z[5]);
	z[49] = exp(z[51]);
	z[11] = exp(z[11]);
	z[50] = exp(z[52]);
	z[12] = exp(z[12]);
	z[43] = exp(z[43]);
	z[13] = exp(z[13]);
	z[44] = exp(z[44]);
	z[15] = exp(z[15]);
	z[45] = exp(z[45]);
	z[7] = exp(z[7]);
	z[46] = exp(z[46]);
	z[8] = exp(z[8]);
	z[47] = exp(z[47]);
	z[9] = exp(z[9]);
	z[48] = exp(z[48]);
	z[6] = exp(z[6]);
	z[16] = h + z[16];
	z[17] = h + z[17];
	z[18] = h + z[18];
	z[19] = h + z[19];
	z[20] = h + z[20];
	z[21] = h + z[21];
	z[22] = h + z[22];
	z[23] = h + z[23];
	z[14] = h + z[14];
	z[26] = h + z[26];
	z[25] = h + z[25];
	z[30] = h + z[30];
	z[31] = h + z[31];
	z[33] = h + z[33];
	z[24] = h + z[24];
	z[34] = h + z[34];
	z[32] = h + z[32];
	z[35] = h + z[35];
	z[36] = h + z[36];
	z[37] = h + z[37];
	z[38] = h + z[38];
	z[39] = h + z[39];
	z[27] = h + z[27];
	z[40] = h + z[40];
	z[28] = h + z[28];
	z[29] = h + z[29];
	z[41] = A*z[41];
	z[42] = A*z[42];
	z[49] = A*z[49];
	z[50] = A*z[50];
	z[43] = A*z[43];
	z[44] = A*z[44];
	z[45] = A*z[45];
	z[46] = A*z[46];
	z[47] = A*z[47];
	z[48] = A*z[48];
	z[16] = z[16]*z[16];
	z[17] = z[17]*z[17];
	z[18] = z[18]*z[18];
	z[19] = z[19]*z[19];
	z[20] = z[20]*z[20];
	z[21] = z[21]*z[21];
	z[22] = z[22]*z[22];
	z[23] = z[23]*z[23];
	z[14] = z[14]*z[14];
	z[26] = z[26]*z[26];
	z[25] = z[25]*z[25];
	z[30] = z[30]*z[30];
	z[31] = z[31]*z[31];
	z[33] = z[33]*z[33];
	z[24] = z[24]*z[24];
	z[34] = z[34]*z[34];
	z[32] = z[32]*z[32];
	z[35] = z[35]*z[35];
	z[36] = z[36]*z[36];
	z[37] = z[37]*z[37];
	z[38] = z[38]*z[38];
	z[39] = z[39]*z[39];
	z[27] = z[27]*z[27];
	z[40] = z[40]*z[40];
	z[28] = z[28]*z[28];
	z[29] = z[29]*z[29];
	z[16] = z[16] + z[3];
	z[17] = z[17] + z[3];
	z[18] = z[18] + z[3];
	z[19] = z[19] + z[3];
	z[20] = z[20] + z[3];
	z[21] = z[21] + z[3];
	z[22] = z[22] + z[3];
	z[23] = z[23] + z[3];
	z[14] = z[14] + z[3];
	z[26] = z[26] + z[3];
	z[25] = z[25] + z[3];
	z[30] = z[3] + z[30];
	z[31] = z[3] + z[31];
	z[33] = z[3] + z[33];
	z[24] = z[24] + z[3];
	z[34] = z[3] + z[34];
	z[32] = z[3] + z[32];
	z[35] = z[3] + z[35];
	z[36] = z[3] + z[36];
	z[37] = z[3] + z[37];
	z[38] = z[3] + z[38];
	z[39] = z[3] + z[39];
	z[27] = z[27] + z[3];
	z[40] = z[3] + z[40];
	z[28] = z[28] + z[3];
	z[3] = z[29] + z[3];
	z[16] = 1./z[16];
	z[17] = 1./z[17];
	z[18] = 1./z[18];
	z[19] = 1./z[19];
	z[20] = 1./z[20];
	z[21] = 1./z[21];
	z[22] = 1./z[22];
	z[23] = 1./z[23];
	z[14] = 1./z[14];
	z[26] = 1./z[26];
	z[25] = 1./z[25];
	z[29] = 1./z[30];
	z[30] = 1./z[31];
	z[31] = 1./z[33];
	z[24] = 1./z[24];
	z[33] = 1./z[34];
	z[32] = 1./z[32];
	z[34] = 1./z[35];
	z[35] = 1./z[36];
	z[36] = 1./z[37];
	z[37] = 1./z[38];
	z[38] = 1./z[39];
	z[27] = 1./z[27];
	z[39] = 1./z[40];
	z[28] = 1./z[28];
	z[3] = 1./z[3];
	z[1] = -z[1];
	z[16] = z[1]*z[16];
	z[17] = z[1]*z[17];
	z[18] = z[1]*z[18];
	z[19] = z[1]*z[19];
	z[20] = z[1]*z[20];
	z[21] = z[1]*z[21];
	z[22] = z[1]*z[22];
	z[23] = z[1]*z[23];
	z[14] = z[1]*z[14];
	z[26] = z[1]*z[26];
	z[25] = z[1]*z[25];
	z[29] = z[1]*z[29];
	z[30] = z[1]*z[30];
	z[31] = z[1]*z[31];
	z[24] = z[1]*z[24];
	z[33] = z[1]*z[33];
	z[32] = z[1]*z[32];
	z[34] = z[1]*z[34];
	z[35] = z[1]*z[35];
	z[36] = z[1]*z[36];
	z[37] = z[1]*z[37];
	z[38] = z[1]*z[38];
	z[27] = z[1]*z[27];
	z[39] = z[1]*z[39];
	z[28] = z[1]*z[28];
	z[1] = z[1]*z[3];
	z[2] = 1. + z[2];
	z[3] = z[16] + z[2];
	z[16] = z[17] + z[2];
	z[17] = z[18] + z[2];
	z[18] = z[19] + z[2];
	z[19] = z[2] + z[20];
	z[20] = z[2] + z[21];
	z[21] = z[2] + z[22];
	z[22] = z[2] + z[23];
	z[14] = z[14] + z[2];
	z[23] = z[2] + z[26];
	z[25] = z[2] + z[25];
	z[26] = z[2] + z[29];
	z[29] = z[2] + z[30];
	z[30] = z[2] + z[31];
	z[24] = z[2] + z[24];
	z[31] = z[2] + z[33];
	z[32] = z[2] + z[32];
	z[33] = z[2] + z[34];
	z[34] = z[2] + z[35];
	z[35] = z[2] + z[36];
	z[36] = z[2] + z[37];
	z[37] = z[2] + z[38];
	z[27] = z[2] + z[27];
	z[38] = z[2] + z[39];
	z[28] = z[2] + z[28];
	z[1] = z[1] + z[2];
	z[2] = beta*z[3];
	z[3] = beta*z[16];
	z[16] = beta*z[17];
	z[17] = beta*z[18];
	z[18] = beta*z[19];
	z[19] = beta*z[20];
	z[20] = beta*z[21];
	z[21] = beta*z[22];
	z[14] = beta*z[14];
	z[22] = beta*z[23];
	z[23] = beta*z[25];
	z[25] = beta*z[26];
	z[26] = beta*z[29];
	z[29] = beta*z[30];
	z[24] = beta*z[24];
	z[30] = beta*z[31];
	z[31] = beta*z[32];
	z[32] = beta*z[33];
	z[33] = beta*z[34];
	z[34] = beta*z[35];
	z[35] = beta*z[36];
	z[36] = beta*z[37];
	z[27] = beta*z[27];
	z[37] = beta*z[38];
	z[28] = beta*z[28];
	z[1] = beta*z[1];
	z[2] = pow(z[2],n);
	z[3] = pow(z[3],n);
	z[17] = z[17] + z[19] + z[22];
	z[16] = z[16] + z[21] + z[25];
	z[14] = z[14] + z[23] + z[26];
	z[18] = z[18] + z[20] + z[29];
	z[1] = z[1] + z[31] + z[33];
	z[19] = z[24] + z[32] + z[36];
	z[20] = z[27] + z[28] + z[35];
	z[21] = z[30] + z[34] + z[37];
	z[2] = 1. + z[2];
	z[3] = 1. + z[3];
	z[17] = pow(z[17],n);
	z[16] = pow(z[16],n);
	z[14] = pow(z[14],n);
	z[18] = pow(z[18],n);
	z[1] = pow(z[1],n);
	z[19] = pow(z[19],n);
	z[20] = pow(z[20],n);
	z[21] = pow(z[21],n);
	z[2] = pow(z[2],z[4]);
	z[3] = pow(z[3],z[4]);
	z[17] = 1. + z[17];
	z[16] = 1. + z[16];
	z[14] = 1. + z[14];
	z[18] = 1. + z[18];
	z[1] = 1. + z[1];
	z[19] = 1. + z[19];
	z[20] = 1. + z[20];
	z[21] = 1. + z[21];
	z[22] = -B*chi;
	z[17] = pow(z[17],z[4]);
	z[16] = pow(z[16],z[4]);
	z[14] = pow(z[14],z[4]);
	z[18] = pow(z[18],z[4]);
	z[1] = pow(z[1],z[4]);
	z[19] = pow(z[19],z[4]);
	z[20] = pow(z[20],z[4]);
	z[4] = pow(z[21],z[4]);
	z[2] = z[10]*z[2]*z[22];
	z[3] = z[22]*z[3]*z[6];
	z[5] = z[17]*z[22]*z[5];
	z[6] = z[11]*z[16]*z[22];
	z[10] = z[12]*z[18]*z[22];
	z[11] = z[13]*z[14]*z[22];
	z[12] = z[15]*z[19]*z[22];
	z[4] = z[22]*z[4]*z[7];
	z[1] = z[1]*z[22]*z[8];
	z[7] = z[20]*z[22]*z[9];
	z[2] = z[2] + z[3] + z[41] + z[48];
	z[2] = 0.5*z[2];
	z[1] = z[1] + z[10] + z[11] + z[12] + z[4] + z[42] + z[49] + z[5] + z[6] + z[7];
	z[1] = z[1] + z[43] + z[44] + z[45] + z[46] + z[47] + z[50];
	z[1] = 0.25*z[1];
	z[1] = z[1] + z[2];

	/* output */
	return z[1];
}