/**
 * @file VarsGenerator.cpp
 * Automatically generate variables for paralleltope generator functions.
 * For high dimensions declaring manually the variables can be tedious...
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "VarsGenerator.h"

/**
 * Constructor that instantiates the variable generator
 *
 * @param[in] dim dimension of the model/parallelotope
 */
VarsGenerator::VarsGenerator(int dim){

	this->dim = dim;

	symbol q1("q1"), q2("q2"), q3("q3"), q4("q4"), q5("q5"), q6("q6"), q7("q7"), q8("q8"), q9("q9"), q10("q10"), q11("q11"), q12("q12"), q13("q13"), q14("q14"), q15("q15"), q16("q16"), q17("q17"), q18("q18"), q19("q19"), q20("q20");
	symbol a1("a1"), a2("a2"), a3("a3"), a4("a4"), a5("a5"), a6("a6"), a7("a7"), a8("a8"), a9("a9"), a10("a10"), a11("a11"), a12("a12"), a13("a13"), a14("a14"), a15("a15"), a16("a16"), a17("a17"), a18("a18"), a19("a19"), a20("a20");
	symbol b1("b1"), b2("b2"), b3("b3"), b4("b4"), b5("b5"), b6("b6"), b7("b7"), b8("b8"), b9("b9"), b10("b10"), b11("b11"), b12("b12"), b13("b13"), b14("b14"), b15("b15"), b16("b16"), b17("b17"), b18("b18"), b19("b19"), b20("b20");
	symbol l1("l1"), l2("l2"), l3("l3"), l4("l4"), l5("l5"), l6("l6"), l7("l7"), l8("l8"), l9("l9"), l10("l10"), l11("l11"), l12("l12"), l13("l13"), l14("l14"), l15("l15"), l16("l16"), l17("l17"), l18("l18"), l19("l19"), l20("l20");

	lst qs, as, bs, ls;
	qs = {q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18,q19,q20};
	as = {a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20};
	bs = {b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20};
	ls = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20};

	symbol u011("u011"), u012("u012"), u013("u013"), u014("u014"), u015("u015"), u016("u016"), u017("u017"), u018("u018"), u019("u019"), u0110("u0110"), u0111("u0111"), u0112("u0112"), u0113("u0113"), u0114("u0114"), u0115("u0115"), u0116("u0116"), u0117("u0117"), u0118("u0118"), u0119("u0119"), u0120("u0120");
	symbol u21("u21"), u22("u22"), u23("u23"), u24("u24"), u25("u25"), u26("u26"), u27("u27"), u28("u28"), u29("u29"), u210("u210"), u211("u211"), u212("u212"), u213("u213"), u214("u214"), u215("u215"), u216("u216"), u217("u217"), u218("u218"), u219("u219"), u220("u220");
	symbol u31("u31"), u32("u32"), u33("u33"), u34("u34"), u35("u35"), u36("u36"), u37("u37"), u38("u38"), u39("u39"), u310("u310"), u311("u311"), u312("u312"), u313("u313"), u314("u314"), u315("u315"), u316("u316"), u317("u317"), u318("u318"), u319("u319"), u320("u320");
	symbol u41("u41"), u42("u42"), u43("u43"), u44("u44"), u45("u45"), u46("u46"), u47("u47"), u48("u48"), u49("u49"), u410("u410"), u411("u411"), u412("u412"), u413("u413"), u414("u414"), u415("u415"), u416("u416"), u417("u417"), u418("u418"), u419("u419"), u420("u420");
	symbol u51("u51"), u52("u52"), u53("u53"), u54("u54"), u55("u55"), u56("u56"), u57("u57"), u58("u58"), u59("u59"), u510("u510"), u511("u511"), u512("u512"), u513("u513"), u514("u514"), u515("u515"), u516("u516"), u517("u517"), u518("u518"), u519("u519"), u520("u520");
	symbol u61("u61"), u62("u62"), u63("u63"), u64("u64"), u65("u65"), u66("u66"), u67("u67"), u68("u68"), u69("u69"), u610("u610"), u611("u611"), u612("u612"), u613("u613"), u614("u614"), u615("u615"), u616("u616"), u617("u617"), u618("u618"), u619("u619"), u620("u620");
	symbol u71("u71"), u72("u72"), u73("u73"), u74("u74"), u75("u75"), u76("u76"), u77("u77"), u78("u78"), u79("u79"), u710("u710"), u711("u711"), u712("u712"), u713("u713"), u714("u714"), u715("u715"), u716("u716"), u717("u717"), u718("u718"), u719("u719"), u720("u720");
	symbol u81("u81"), u82("u82"), u83("u83"), u84("u84"), u85("u85"), u86("u86"), u87("u87"), u88("u88"), u89("u89"), u810("u810"), u811("u811"), u812("u812"), u813("u813"), u814("u814"), u815("u815"), u816("u816"), u817("u817"), u818("u818"), u819("u819"), u820("u820");
	symbol u91("u91"), u92("u92"), u93("u93"), u94("u94"), u95("u95"), u96("u96"), u97("u97"), u98("u98"), u99("u99"), u910("u910"), u911("u911"), u912("u912"), u913("u913"), u914("u914"), u915("u915"), u916("u916"), u917("u917"), u918("u918"), u919("u919"), u920("u920");
	symbol u101("u101"), u102("u102"), u103("u103"), u104("u104"), u105("u105"), u106("u106"), u107("u107"), u108("u108"), u109("u109"), u1010("u1010"), u1011("u1011"), u1012("u1012"), u1013("u1013"), u1014("u1014"), u1015("u1015"), u1016("u1016"), u1017("u1017"), u1018("u1018"), u1019("u1019"), u1020("u1020");
	symbol u111("u111"), u112("u112"), u113("u113"), u114("u114"), u115("u115"), u116("u116"), u117("u117"), u118("u118"), u119("u119"), u1110("u1110"), u1111("u1111"), u1112("u1112"), u1113("u1113"), u1114("u1114"), u1115("u1115"), u1116("u1116"), u1117("u1117"), u1118("u1118"), u1119("u1119"), u1120("u1120");
	symbol u121("u121"), u122("u122"), u123("u123"), u124("u124"), u125("u125"), u126("u126"), u127("u127"), u128("u128"), u129("u129"), u1210("u1210"), u1211("u1211"), u1212("u1212"), u1213("u1213"), u1214("u1214"), u1215("u1215"), u1216("u1216"), u1217("u1217"), u1218("u1218"), u1219("u1219"), u1220("u1220");
	symbol u131("u131"), u132("u132"), u133("u133"), u134("u134"), u135("u135"), u136("u136"), u137("u137"), u138("u138"), u139("u139"), u1310("u1310"), u1311("u1311"), u1312("u1312"), u1313("u1313"), u1314("u1314"), u1315("u1315"), u1316("u1316"), u1317("u1317"), u1318("u1318"), u1319("u1319"), u1320("u1320");
	symbol u141("u141"), u142("u142"), u143("u143"), u144("u144"), u145("u145"), u146("u146"), u147("u147"), u148("u148"), u149("u149"), u1410("u1410"), u1411("u1411"), u1412("u1412"), u1413("u1413"), u1414("u1414"), u1415("u1415"), u1416("u1416"), u1417("u1417"), u1418("u1418"), u1419("u1419"), u1420("u1420");
	symbol u151("u151"), u152("u152"), u153("u153"), u154("u154"), u155("u155"), u156("u156"), u157("u157"), u158("u158"), u159("u159"), u1510("u1510"), u1511("u1511"), u1512("u1512"), u1513("u1513"), u1514("u1514"), u1515("u1515"), u1516("u1516"), u1517("u1517"), u1518("u1518"), u1519("u1519"), u1520("u1520");
	symbol u161("u161"), u162("u162"), u163("u163"), u164("u164"), u165("u165"), u166("u166"), u167("u167"), u168("u168"), u169("u169"), u1610("u1610"), u1611("u1611"), u1612("u1612"), u1613("u1613"), u1614("u1614"), u1615("u1615"), u1616("u1616"), u1617("u1617"), u1618("u1618"), u1619("u1619"), u1620("u1620");
	symbol u171("u171"), u172("u172"), u173("u173"), u174("u174"), u175("u175"), u176("u176"), u177("u177"), u178("u178"), u179("u179"), u1710("u1710"), u1711("u1711"), u1712("u1712"), u1713("u1713"), u1714("u1714"), u1715("u1715"), u1716("u1716"), u1717("u1717"), u1718("u1718"), u1719("u1719"), u1720("u1720");

	lst u1s, u2s, u3s, u4s, u5s, u6s, u7s, u8s, u9s, u10s, u11s, u12s, u13s, u14s, u15s, u16s, u17s, u18s, u19s, u20s;
	u1s = {u011, u012, u013, u014, u015, u016, u017, u018, u019, u0110, u0111, u0112, u0113, u0114, u0115, u0116, u0117, u0118, u0119, u0120};
	u2s = {u21, u22, u23, u24, u25, u26, u27, u28, u29, u210, u211, u212, u213, u214, u215, u216, u217, u218, u219, u220};
	u3s = {u31, u32, u33, u34, u35, u36, u37, u38, u39, u310, u311, u312, u313, u314, u315, u316, u317, u318, u319, u320};
	u4s = {u41, u42, u43, u44, u45, u46, u47, u48, u49, u410, u411, u412, u413, u414, u415, u416, u417, u418, u419, u420};
	u5s = {u51, u52, u53, u54, u55, u56, u57, u58, u59, u510, u511, u512, u513, u514, u515, u516, u517, u518, u519, u520};
	u6s = {u61, u62, u63, u64, u65, u66, u67, u68, u69, u610, u611, u612, u613, u614, u615, u616, u617, u618, u619, u620};
	u7s = {u71, u72, u73, u74, u75, u76, u77, u78, u79, u710, u711, u712, u713, u714, u715, u716, u717, u718, u719, u720};
	u8s = {u81, u82, u83, u84, u85, u86, u87, u88, u89, u810, u811, u812, u813, u814, u815, u816, u817, u818, u819, u820};
	u9s = {u91, u92, u93, u94, u95, u96, u97, u98, u99, u910, u911, u912, u913, u914, u915, u916, u917, u918, u919, u920};
	u10s = {u101, u102, u103, u104, u105, u106, u107, u108, u109, u1010, u1011, u1012, u1013, u1014, u1015, u1016, u1017, u1018, u1019, u1020};
	u11s = {u111, u112, u113, u114, u115, u116, u117, u118, u119, u1110, u1111, u1112, u1113, u1114, u1115, u1116, u1117, u1118, u1119, u1120};
	u12s = {u121, u122, u123, u124, u125, u126, u127, u128, u129, u1210, u1211, u1212, u1213, u1214, u1215, u1216, u1217, u1218, u1219, u1220};
	u13s = {u131, u132, u133, u134, u135, u136, u137, u138, u139, u1310, u1311, u1312, u1313, u1314, u1315, u1316, u1317, u1318, u1319, u1320};
	u14s = {u141, u142, u143, u144, u145, u146, u147, u148, u149, u1410, u1411, u1412, u1413, u1414, u1415, u1416, u1417, u1418, u1419, u1420};
	u15s = {u151, u152, u153, u154, u155, u156, u157, u158, u159, u1510, u1511, u1512, u1513, u1514, u1515, u1516, u1517, u1518, u1519, u1520};
	u16s = {u161, u162, u163, u164, u165, u166, u167, u168, u169, u1610, u1611, u1612, u1613, u1614, u1615, u1616, u1617, u1618, u1619, u1620};
	u17s = {u171, u172, u173, u174, u175, u176, u177, u178, u179, u1710, u1711, u1712, u1713, u1714, u1715, u1716, u1717, u1718, u1719, u1720};

	vector< lst > us;
	us.push_back(u1s); us.push_back(u2s); us.push_back(u3s);
	us.push_back(u4s); us.push_back(u5s); us.push_back(u6s);
	us.push_back(u7s); us.push_back(u8s); us.push_back(u9s);
	us.push_back(u10s); us.push_back(u11s); us.push_back(u12s);
	us.push_back(u13s); us.push_back(u14s); us.push_back(u15s);
	us.push_back(u16s); us.push_back(u17s);


	for(int i=0; i<this->dim; i++){
		this->qs.append(qs[i]);
		this->as.append(as[i]);
		this->bs.append(bs[i]);
		this->ls.append(ls[i]);
	}

	for(int i=0; i<this->dim; i++){
		lst us_i;
		for(int j=0; j<this->dim; j++){
			us_i.append(us[i][j]);
		}
		this->us.push_back(us_i);
	}

}

/**
 * Generate a box out of the variables
 *
 * @param[in] b box offsets
 * @returns generatred box
 */
LinearSystem* VarsGenerator::genBox(vector<double> b){
	int n = b.size()/2;
	vector<double> Ai(n,0);
	vector< vector<double> > A(2*n,Ai);

	for(int i=0; i<n; i++){
		A[2*i][i] = 1;
		A[2*i+1][i] = -1;
	}

	LinearSystem *Ab = new LinearSystem(A,b);
	return Ab;
}


VarsGenerator::~VarsGenerator() {
	// TODO Auto-generated destructor stub
}
