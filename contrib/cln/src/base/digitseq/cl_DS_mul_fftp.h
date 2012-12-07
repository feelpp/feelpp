// Fast integer multiplication using FFT in a modular ring.
// Bruno Haible 5.5.1996

// FFT in the complex domain has the drawback that it needs careful round-off
// error analysis. So here we choose another field of characteristic 0: Q_p.
// Since Q_p contains exactly the (p-1)th roots of unity, we choose
// p == 1 mod N and have the Nth roots of unity (N = 2^n) in Q_p and
// even in Z_p. Actually, we compute in Z/(p^m Z).

// All operations the FFT algorithm needs is addition, subtraction,
// multiplication, multiplication by the Nth root of unity and division
// by N. Hence we can use the domain Z/(p^m Z) even if p is not a prime!

// We want to compute the convolution of N 32-bit words. The resulting
// words are < (2^32)^2 * N. If is safe to compute in Z/pZ with p = 2^94 + 1
// or p = 7*2^92 + 1. We choose p < 2^95 so that we can easily represent every
// element of Z/pZ as three 32-bit words.

#if !(intDsize==32)
#error "fft mod p implemented only for intDsize==32"
#endif

#if 0
typedef union {
	#if CL_DS_BIG_ENDIAN_P
	struct { uint32 w2; uint32 w1; uint32 w0; };
	#else
	struct { uint32 w0; uint32 w1; uint32 w2; };
	#endif
	uintD _w[3];
} fftp_word;
#else
// typedef struct { uint32 w2; uint32 w1; uint32 w0; } fftp_word;
// typedef struct { uint32 w0; uint32 w1; uint32 w2; } fftp_word;
typedef struct { uintD _w[3]; } fftp_word;
#endif
#if CL_DS_BIG_ENDIAN_P
  #define w2 _w[0]
  #define w1 _w[1]
  #define w0 _w[2]
  #define W3(W2,W1,W0)  { W2, W1, W0 }
#else
  #define w0 _w[0]
  #define w1 _w[1]
  #define w2 _w[2]
  #define W3(W2,W1,W0)  { W0, W1, W2 }
#endif

#if 0
// p = 19807040628566084398385987585 = 5 * 3761 * 7484047069 * 140737471578113
static const fftp_word p = W3( 1L<<30, 0, 1 ); // p = 2^94 + 1
#define FFT_P_94
static const fftp_word fftp_roots_of_1 [24+1] =
  // roots_of_1[n] is a (2^n)th root of unity in Z/pZ.
  // (Also roots_of_1[n-1] = roots_of_1[n]^2, but we don't need this.)
  // (To build this table, you need to compute roots of unity modulo the
  // factors of p and combine them using the Chinese Remainder Theorem.
  // Or ask me for "quadmod.lsp".)
  {
    W3( 0x00000000, 0x00000000, 0x00000001 ), //                             1
    W3( 0x0000003F, 0xFFFFFFFF, 0xFF800000 ), //        1180591620717402914816
    W3( 0x20000040, 0x00004000, 0x00000001 ), //  9903521494874733285348474881
    W3( 0x3688E9A7, 0xDD78E2A9, 0x1E75974D ), // 16877707849775746711303853901
    W3( 0x286E6589, 0x5E86C1E0, 0x42710379 ), // 12512861726041464545960067961
    W3( 0x00D79325, 0x1A884885, 0xEA46D6C5 ), //   260613923531515619478787781
    W3( 0x1950B480, 0xC387CEE5, 0xA69C443F ), //  7834691712342412468047070271
    W3( 0x19DC9D08, 0x11CADC6A, 0x5BA8B123 ), //  8003830486242687653832601891
    W3( 0x21D6D905, 0xB8BAC7C3, 0xC3841613 ), // 10472740308573592285123712531
    W3( 0x27D73986, 0x6AF6BD27, 0x7A6D7909 ), // 12330106088710388189231937801
    W3( 0x20D4698B, 0x0039D457, 0xA092AECF ), // 10160311000635748689099534031
    W3( 0x049BD1C4, 0xA94F001A, 0xFA76E358 ), //  1426314143682376031341568856
    W3( 0x26DD7228, 0x09400257, 0x9BB49CB9 ), // 12028142067661291067236719801
    W3( 0x12DAA9AD, 0xAF9435A9, 0xD50FF483 ), //  5835077289334326375656453251
    W3( 0x0B7CDA03, 0x9418702E, 0x7CD934CA ), //  3555271451571910239441204426
    W3( 0x2D272FCF, 0xB8644522, 0x68EAD40B ), // 13974199331913037576372147211
    W3( 0x00EDA06E, 0x0114DA26, 0xE8D84BA9 ), //   287273027105701319912475561
    W3( 0x2219C2C4, 0xFD3145C6, 0xDD019359 ), // 10553633252320053510122083161
    W3( 0x1764F007, 0x4F5D5FD4, 0xDAB10AFC ), //  7240181310654329198595869436
    W3( 0x01AA13EE, 0x2D1CD906, 0x11D5B1EB ), //   515096517694807745704079851
    W3( 0x27038944, 0x5A37BAAD, 0x5CECA64C ), // 12074190385578921562318087756
    W3( 0x2459CF22, 0xF625FD38, 0xADB48511 ), // 11250032926302238120667809041
    W3( 0x25B6C6A8, 0xD684063F, 0x7ABAD1EF ), // 11671908005633729316324561391
    W3( 0x1C1A2BC6, 0x12B253F1, 0x0D1BBCB7 ), //  8697219061868963805380983991
    W3( 0x198F3FE2, 0x5EE9919F, 0x535E80D5 }  //  7910303322630257758732976341
  };
// Sadly, this p doesn't work because we don't find a (2^n)th root of unity w
// such that  w^(2^(n-1)) = -1  mod p. However, our algorithm below assumes
// that  w^(2^(n-1)) = -1...
#else
// p = 34662321099990647697175478273, a prime
static const fftp_word p = W3( 7L<<28, 0, 1 ); // p = 7 * 2^92 + 1
#define FFT_P_92
static const fftp_word fftp_roots_of_1 [92+1] =
  // roots_of_1[n] is a (2^n)th root of unity in Z/pZ.
  // (Also roots_of_1[n-1] = roots_of_1[n]^2, but we don't need this.)
  {
    W3( 0x00000000, 0x00000000, 0x00000001 ), //                             1
    W3( 0x70000000, 0x00000000, 0x00000000 ), // 34662321099990647697175478272
    W3( 0x064AF70F, 0x997E62CE, 0x77953100 ), //  1947537281862369253065568512
    W3( 0x261B8E96, 0xC3AD4296, 0xDA1BFA93 ), // 11793744727492885369350519443
    W3( 0x096EA949, 0x6EDCAF05, 0x47C92A4F ), //  2919146363086089454841571919
    W3( 0x366A8C3F, 0x7BF1436D, 0x2333BE9E ), // 16840998969615256469762195102
    W3( 0x27569FA8, 0xAE1775F1, 0xB21956A0 ), // 12174636971387721414084220576
    W3( 0x16CABB8B, 0xBAA59813, 0x62FCBCD9 ), //  7053758891710792545762852057
    W3( 0x1AE130A3, 0xF909B101, 0xB6BA30CF ), //  8318848263123793919933558991
    W3( 0x32AE8FEE, 0x6B1A656B, 0xED02BF24 ), // 15685283280129931240441823012
    W3( 0x2D1EE047, 0x5AEDC882, 0x8E96BCCC ), // 13964152342912072497719852236
    W3( 0x222A18FD, 0x3BF40635, 0xBFDEA8AD ), // 10573383226491471052459124909
    W3( 0x10534EE6, 0xED5A55D4, 0x06AE2155 ), //  5052473604609413010647032149
    W3( 0x02F3BFA3, 0x2D816786, 0xE6C27B3C ), //   913643975905572976593107772
    W3( 0x0B0CD0A5, 0x9A1FF4F7, 0x2624A5E1 ), //  3419827524917244902802499041
    W3( 0x257A492F, 0x156C141C, 0xFC5D75F4 ), // 11598779914676604137587439092
    W3( 0x061FB92A, 0xB1A1F41A, 0x7006920F ), //  1895261184698485907279745551
    W3( 0x2A4E1471, 0xDDB96073, 0xD8DDBB71 ), // 13092763174215078887900953457
    W3( 0x213B469E, 0xD72A84CA, 0xAAA477F2 ), // 10284665443205365583657072626
    W3( 0x1D7EF67C, 0x3DC2DA37, 0x4C86E9DC ), //  9128553932091860654576036316
    W3( 0x0CB7AA67, 0x2E087ED8, 0x2675D6E3 ), //  3935858248479385987820410595
    W3( 0x00BD7B24, 0x68388052, 0x57FFFB10 ), //   229068502577238003716979472
    W3( 0x1E1724A6, 0xBA587C3D, 0x0C12825B ), //  9312528669272006966417457755
    W3( 0x20595EF0, 0xC89DA33B, 0x3CB5583B ), // 10011563056352601486430394427
    W3( 0x15E730B2, 0x6D34E9EB, 0x71CCE555 ), //  6778677035560020292206912853
    W3( 0x015EFDBB, 0xC0A80C3B, 0xE4B1E017 ), //   424322259008787317821399063
    W3( 0x1B81FC63, 0x0C694944, 0x8EB481BF ), //  8513238559382277026756198847
    W3( 0x1AF53421, 0x5DCAA1A4, 0xD0C15A03 ), //  8343043259718611508685527555
    W3( 0x2F2B6B58, 0xBB60E464, 0x37A7DE2E ), // 14598286201875835624993840686
    W3( 0x27B4AB13, 0x54617640, 0xE86E757A ), // 12288329911800070034603013498
    W3( 0x041A31D2, 0xF0AC8E3C, 0x8AA4FD27 ), //  1269607397711669380834983207
    W3( 0x1A52F484, 0x39AC5917, 0x34E3F1F7 ), //  8146896869111203814625767927
    W3( 0x048FC120, 0x50F6ECBF, 0x268D86A8 ), //  1411728444351387120148776616
    W3( 0x27A2C427, 0x001F1239, 0x93380047 ), // 12266687669072434694473646151
    W3( 0x2E7E8DFB, 0x2411A754, 0xE12A9B1D ), // 14389305591459206001391737629
    W3( 0x29F14702, 0x40B3E1E2, 0xF7D71A8D ), // 12980571854778363745245010573
    W3( 0x3158DCE7, 0x8B8FEB32, 0x1DE35D24 ), // 15272194145252623177165790500
    W3( 0x12484C07, 0x437ED373, 0x9E45F602 ), //  5658131869639928287764805122
    W3( 0x1AEAE06E, 0xB905C908, 0x4389BF5F ), //  8330558749711089231534341983
    W3( 0x27BC0045, 0x43024FEB, 0xEC880258 ), // 12297194714773858676269122136
    W3( 0x2EFE1CBC, 0x0D2FAA94, 0xB4EA69A6 ), // 14543513305163560781242591654
    W3( 0x0B0D3D8B, 0xD779F105, 0x920367FA ), //  3420341787669373425792804858
    W3( 0x2D4D7BA9, 0x0970D8CF, 0x8CE6D7EC ), // 14020496699328277892009744364
    W3( 0x00DC5971, 0x0209470E, 0x713F2B27 ), //   266386055561000736260041511
    W3( 0x27E54E26, 0x53BA0137, 0xDD6740B3 ), // 12347128447319282384829366451
    W3( 0x2143A889, 0x8F2B57F5, 0xFB8181C1 ), // 10294799249108063706647986625
    W3( 0x1125419F, 0x5C4E0608, 0xE0AC0396 ), //  5306285315793562029414679446
    W3( 0x15B61D90, 0x63A27BB0, 0x26402B32 ), //  6719349317556695539371748146
    W3( 0x03B582FC, 0x419EF656, 0xB06BBC35 ), //  1147889163765050226454019125
    W3( 0x08FF62E1, 0xA3BB1145, 0xDA998F77 ), //  2784623116803271439773437815
    W3( 0x101978AF, 0xF93CBFA1, 0xB788B5A3 ), //  4982553232749484200897852835
    W3( 0x061334DE, 0x8FE5C6E9, 0x2B2309D6 ), //  1880129318103954373583505878
    W3( 0x343C6E7C, 0x8019BB43, 0xD954E744 ), // 16166277816826816936484857668
    W3( 0x06506A03, 0x0E6DE333, 0xF8011494 ), //  1954124751724394051182597268
    W3( 0x34892A42, 0x6502DAA3, 0x8FDA6971 ), // 16259042912153157504364865905
    W3( 0x0EF2C4BD, 0xF42D9711, 0xC32CEA49 ), //  4626279273705729025744104009
    W3( 0x24511305, 0x4F1EAE2C, 0x62FB10F4 ), // 11239473167855288013010178292
    W3( 0x14E5A052, 0xF1748A9C, 0xDD536730 ), //  6467301317787608309692589872
    W3( 0x0621D0A7, 0x0A5188AF, 0x7316C352 ), //  1897789944553576071437927250
    W3( 0x234498F0, 0xDF078E95, 0x6FEED50B ), // 10914904542475816633386325259
    W3( 0x029E4925, 0x948D6D57, 0xD4DF93A6 ), //   810325725128913871737688998
    W3( 0x11BB3805, 0x0589D746, 0x852F3E2F ), //  5487578840386649632552205871
    W3( 0x1D4370CA, 0xA4441B85, 0xC9606FE0 ), //  9056595957858187419376971744
    W3( 0x1C536F7D, 0x77D44926, 0x8DDB8932 ), //  8766447615182890705620797746
    W3( 0x3498CE71, 0xB726A4D3, 0xF4F3C813 ), // 16277952140466335672647796755
    W3( 0x1E4A297E, 0xAC13196E, 0xFACD8102 ), //  9374206759006667727054930178
    W3( 0x0E7C2CCC, 0xC940C98B, 0x0BC0CA49 ), //  4482908500893894680116251209
    W3( 0x124CF912, 0xD84438FD, 0x9C03585F ), //  5663784755954194195257972831
    W3( 0x06180FF8, 0xD447BEBE, 0xDB8821E7 ), //  1885999704184999223512015335
    W3( 0x1ED2EB11, 0x0687EC7C, 0xBE3436C8 ), //  9539534786948152514714023624
    W3( 0x30EBB35C, 0x59616A3C, 0x502CBB52 ), // 15140225046175435352009653074
    W3( 0x33E24883, 0xEDA36D60, 0xA25C8E5F ), // 16057295180155395438855097951
    W3( 0x0D879ED9, 0x076BAB06, 0x9BE12AA2 ), //  4187260250707927256570866338
    W3( 0x1A1B6C9C, 0x0966383B, 0x54123A87 ), //  8079764146434082816365050503
    W3( 0x31BD863A, 0xA2A6505C, 0xD759E6CF ), // 15393886339893077529104869071
    W3( 0x3209AF0A, 0x5E5055A1, 0x480AF03F ), // 15485957428841754012708171839
    W3( 0x1A4CC03C, 0xC8AA650B, 0x7F4DBCE9 ), //  8139396433274519652257348841
    W3( 0x3596471F, 0xB99D2EA3, 0xA3433E0C ), // 16584380266717730797139475980
    W3( 0x28E87642, 0x98E21FCE, 0xDE1B53EA ), // 12660429650750886812340409322
    W3( 0x20161DB9, 0xCDC199E9, 0x0A6BEDF2 ), //  9930257058416521571476434418
    W3( 0x1D0DC095, 0x2C40D22B, 0x088549BA ), //  8991690766592354745340742074
    W3( 0x2FCC953C, 0xA8B62408, 0x50FC4C29 ), // 14793121080372138443684989993
    W3( 0x0F854B39, 0xF659B4B5, 0xD2B0A6AC ), //  4803417528030967235217499820
    W3( 0x30E087D4, 0x02F3BBAB, 0xBA503373 ), // 15126721285415891216108630899
    W3( 0x0DAF660C, 0x26B99C42, 0x98B8BE05 ), //  4235349051642660841298902533
    W3( 0x0ED6AE0E, 0xCD02982A, 0xD233F0D9 ), //  4592322227691334993146278105
    W3( 0x3415EB9B, 0x4B61C19F, 0xB21F1255 ), // 16119720573722492095181034069
    W3( 0x1015A729, 0x20A1FAA2, 0x0D094529 ), //  4977936993224010619482096937
    W3( 0x1D2E3AD2, 0x7093579F, 0x1C93C97B ), //  9030953651705465548198627707
    W3( 0x130EAA8F, 0x859C980F, 0xD9E7E8ED ), //  5897945597894388791627999469
    W3( 0x2B7CA1C8, 0xFC34C5B5, 0x9C0B1C0C ), // 13458526232475976507763399692
    W3( 0x22367055, 0xA53B526A, 0x7505EABE ), // 10588302813110450634719881918
    W3( 0x344FEF55, 0x0B77067F, 0x38999E77 )  // 16189855864848287589134343799
  };
#endif

// Define this if you want the external loops instead of inline operations.
#define FFTP_EXTERNAL_LOOPS

// Define this for (cheap) consistency checks.
//#define DEBUG_FFTP

// Define this for extensive consistency checks.
//#define DEBUG_FFTP_OPERATIONS

// Define the algorithm of the backward FFT:
// Either FORWARD (a normal FFT followed by a permutation)
// or     RECIPROOT (an FFT with reciprocal root of unity)
// or     CLEVER (an FFT with reciprocal root of unity but clever computation
//                of the reciprocals).
// Drawback of FORWARD: the permutation pass.
// Drawback of RECIPROOT: need all the powers of the root, not only half of them.
#define FORWARD   42
#define RECIPROOT 43
#define CLEVER    44
#define FFTP_BACKWARD CLEVER

// r := a + b
static inline void add (const fftp_word& a, const fftp_word& b, fftp_word& r)
{
#ifdef FFTP_EXTERNAL_LOOPS
	add_loop_lsp(arrayLSDptr(a._w,3),arrayLSDptr(b._w,3),arrayLSDptr(r._w,3),3);
#else
	var uint32 tmp;

	tmp = a.w0 + b.w0;
	if (tmp >= a.w0) {
		// no carry
		r.w0 = tmp;
		tmp = a.w1 + b.w1;
		if (tmp >= a.w1) goto no_carry_1; else goto carry_1;
	} else {
		// carry
		r.w0 = tmp;
		tmp = a.w1 + b.w1 + 1;
		if (tmp > a.w1) goto no_carry_1; else goto carry_1;
	}
	if (1) {
		no_carry_1: // no carry
		r.w1 = tmp;
		tmp = a.w2 + b.w2;
	} else {
		carry_1: // carry
		r.w1 = tmp;
		tmp = a.w2 + b.w2 + 1;
	}
	r.w2 = tmp;
#endif
}

// r := a - b
static inline void sub (const fftp_word& a, const fftp_word& b, fftp_word& r)
{
#ifdef FFTP_EXTERNAL_LOOPS
	sub_loop_lsp(arrayLSDptr(a._w,3),arrayLSDptr(b._w,3),arrayLSDptr(r._w,3),3);
#else
	var uint32 tmp;

	tmp = a.w0 - b.w0;
	if (tmp <= a.w0) {
		// no carry
		r.w0 = tmp;
		tmp = a.w1 - b.w1;
		if (tmp <= a.w1) goto no_carry_1; else goto carry_1;
	} else {
		// carry
		r.w0 = tmp;
		tmp = a.w1 - b.w1 - 1;
		if (tmp < a.w1) goto no_carry_1; else goto carry_1;
	}
	if (1) {
		no_carry_1: // no carry
		r.w1 = tmp;
		tmp = a.w2 - b.w2;
	} else {
		carry_1: // carry
		r.w1 = tmp;
		tmp = a.w2 - b.w2 - 1;
	}
	r.w2 = tmp;
#endif
}

// b := a >> 1
static inline void shift (const fftp_word& a, fftp_word& b)
{
#ifdef FFTP_EXTERNAL_LOOPS
	#ifdef DEBUG_FFTP
	if (shiftrightcopy_loop_msp(arrayMSDptr(a._w,3),arrayMSDptr(b._w,3),3,1,0))
		throw runtime_exception();
	#else
	shiftrightcopy_loop_msp(arrayMSDptr(a._w,3),arrayMSDptr(b._w,3),3,1,0);
	#endif
#else
	var uint32 tmp, carry;

	tmp = a.w2;
	b.w2 = a.w2 >> 1;
	carry = tmp << 31;
	tmp = a.w1;
	b.w1 = (tmp >> 1) | carry;
	carry = tmp << 31;
	tmp = a.w0;
	b.w0 = (tmp >> 1) | carry;
	#ifdef DEBUG_FFTP
	carry = tmp << 31;
	if (carry)
		throw runtime_exception();
	#endif
#endif
}

#ifdef DEBUG_FFTP_OPERATIONS
#define check_fftp_word(x)  if (compare_loop_msp(arrayMSDptr((x)._w,3),arrayMSDptr(p._w,3),3) >= 0) throw runtime_exception()
#else
#define check_fftp_word(x)
#endif

// r := (a + b) mod p
static inline void addp (const fftp_word& a, const fftp_word& b, fftp_word& r)
{
	check_fftp_word(a); check_fftp_word(b);
#ifdef FFTP_EXTERNAL_LOOPS
	add(a,b, r);
	if (compare_loop_msp(arrayMSDptr(r._w,3),arrayMSDptr(p._w,3),3) >= 0)
		sub(r,p, r);
#else
	add(a,b, r);
	if ((r.w2 > p.w2)
	    || ((r.w2 == p.w2)
	        && ((r.w1 > p.w1)
	            || ((r.w1 == p.w1)
	                && (r.w0 >= p.w0)))))
		sub(r,p, r);
#endif
	check_fftp_word(r);
}

// r := (a - b) mod p
static inline void subp (const fftp_word& a, const fftp_word& b, fftp_word& r)
{
	check_fftp_word(a); check_fftp_word(b);
	sub(a,b, r);
	if ((sint32)r.w2 < 0)
		add(r,p, r);
	check_fftp_word(r);
}

// r := (a * b) mod p
static void mulp (const fftp_word& a, const fftp_word& b, fftp_word& r)
{
	check_fftp_word(a); check_fftp_word(b);
#if defined(FFT_P_94)
	var uintD c[6];
	var uintD* const cLSDptr = arrayLSDptr(c,6);
	// Multiply the two words, using the standard method.
	mulu_2loop(arrayLSDptr(a._w,3),3, arrayLSDptr(b._w,3),3, cLSDptr);
	// c[0..5] now contains the product.
	// Divide by p.
	// To divide c (0 <= c < p^2) by p = 2^n+1,
	// we set q := floor(c/2^n) and r := c - q*p = (c mod 2^n) - q.
	// If this becomes negative, set r := r + p (at most twice).
	// (This works because  floor(c/p) <= q <= floor(c/p)+2.)
	// (Actually, here, 0 <= c <= (p-1)^2, hence
	// floor(c/p) <= q <= floor(c/p)+1, so we have
	// to set r := r + p at most once!)
	// n = 94 = 3*32-2 = 2*32+30.
	shiftleft_loop_lsp(cLSDptr lspop 3,3,2,lspref(cLSDptr,2)>>30);
	lspref(cLSDptr,2) &= bit(30)-1;
	// c[0..2] now contains q, c[3..5] contains (c mod 2^n).
	#if 0
	if (compare_loop_msp(cLSDptr lspop 6,arrayMSDptr(p._w,3),3) >= 0) // q >= p ?
		subfrom_loop_lsp(arrayLSDptr(p._w,3),cLSDptr lspop 3,3); // q -= p;
	#endif
	if (subfrom_loop_lsp(cLSDptr lspop 3,cLSDptr,3)) // (c mod 2^n) - q
		addto_loop_lsp(arrayLSDptr(p._w,3),cLSDptr,3);
	r.w2 = lspref(cLSDptr,2); r.w1 = lspref(cLSDptr,1); r.w0 = lspref(cLSDptr,0);
#elif defined(FFT_P_92)
	var uintD c[7];
	var uintD* const cLSDptr = arrayLSDptr(c,7);
	// Multiply the two words, using the standard method.
	mulu_2loop(arrayLSDptr(a._w,3),3, arrayLSDptr(b._w,3),3, cLSDptr);
	// c[1..6] now contains the product.
	// Divide by p.
	// To divide c (0 <= c < p^2) by p = 7*2^n+1,
	// we set q := floor(floor(c/2^n)/7) and
	// r := c - q*p = (floor(c/2^n) mod 7)*2^n + (c mod 2^n) - q.
	// If this becomes negative, set r := r + p.
	// (As above, since 0 <= c <= (p-1)^2, we have
	// floor(c/p) <= q <= floor(c/p)+1, so we have
	// to set r := r + p at most once!)
	// n = 92 = 3*32-4 = 2*32+28.
	lspref(cLSDptr,6) = shiftleft_loop_lsp(cLSDptr lspop 3,3,4,lspref(cLSDptr,2)>>28);
	lspref(cLSDptr,2) &= bit(28)-1;
	// c[0..3] now contains floor(c/2^n), c[4..6] contains (c mod 2^n).
	var uintD remainder = divu_loop_msp(7,cLSDptr lspop 7,4);
	lspref(cLSDptr,2) |= remainder << 28;
	// c[0..3] now contains q, c[4..6] contains (c mod 7*2^n).
	#ifdef DEBUG_FFTP
	if (lspref(cLSDptr,6) > 0)
		throw runtime_exception();
	#endif
	#if 0
	if (compare_loop_msp(cLSDptr lspop 6,arrayMSDptr(p._w,3),3) >= 0) // q >= p ?
		subfrom_loop_lsp(arrayLSDptr(p._w,3),cLSDptr lspop 3,3); // q -= p;
	#endif
	if (subfrom_loop_lsp(cLSDptr lspop 3,cLSDptr,3)) // (c mod 2^n) - q
		addto_loop_lsp(arrayLSDptr(p._w,3),cLSDptr,3);
	r.w2 = lspref(cLSDptr,2); r.w1 = lspref(cLSDptr,1); r.w0 = lspref(cLSDptr,0);
#else
#error "mulp not implemented for this prime"
#endif
	if ((sint32)r.w2 < 0)
		throw runtime_exception();
	check_fftp_word(r);
}
#ifdef DEBUG_FFTP_OPERATIONS
static void mulp_doublecheck (const fftp_word& a, const fftp_word& b, fftp_word& r)
{
	fftp_word zero, ma, mb, or;
	subp(a,a, zero);
	subp(zero,a, ma);
	subp(zero,b, mb);
	mulp(ma,mb, or);
	mulp(a,b, r);
	if (compare_loop_msp(arrayMSDptr(r._w,3),arrayMSDptr(or._w,3),3))
		throw runtime_exception();
}
#define mulp mulp_doublecheck
#endif /* DEBUG_FFTP_OPERATIONS */

// b := (a / 2) mod p
static inline void shiftp (const fftp_word& a, fftp_word& b)
{
	check_fftp_word(a);
	if (a.w0 & 1) {
		var fftp_word a_even;
		add(a,p, a_even);
		shift(a_even, b);
	} else
		shift(a, b);
	check_fftp_word(b);
}

#ifndef _BIT_REVERSE
#define _BIT_REVERSE
// Reverse an n-bit number x. n>0.
static uintC bit_reverse (uintL n, uintC x)
{
	var uintC y = 0;
	do {
		y <<= 1;
		y |= (x & 1);
		x >>= 1;
	} while (!(--n == 0));
	return y;
}
#endif

// Compute an convolution mod p using FFT: z[0..N-1] := x[0..N-1] * y[0..N-1].
static void fftp_convolution (const uintL n, const uintC N, // N = 2^n
                              fftp_word * x, // N words
                              fftp_word * y, // N words
                              fftp_word * z  // N words result
                             )
{
	CL_ALLOCA_STACK;
	#if (FFTP_BACKWARD == RECIPROOT) || defined(DEBUG_FFTP)
	var fftp_word* const w = cl_alloc_array(fftp_word,N);
	#else
	var fftp_word* const w = cl_alloc_array(fftp_word,(N>>1)+1);
	#endif
	var uintC i;
	// Initialize w[i] to w^i, w a primitive N-th root of unity.
	w[0] = fftp_roots_of_1[0];
	w[1] = fftp_roots_of_1[n];
	#if (FFTP_BACKWARD == RECIPROOT) || defined(DEBUG_FFTP)
	for (i = 2; i < N; i++)
		mulp(w[i-1],fftp_roots_of_1[n], w[i]);
	#else // need only half of the roots
	for (i = 2; i < N>>1; i++)
		mulp(w[i-1],fftp_roots_of_1[n], w[i]);
	#endif
	#ifdef DEBUG_FFTP
	// Check that w is really a primitive N-th root of unity.
	{
		var fftp_word w_N;
		mulp(w[N-1],fftp_roots_of_1[n], w_N);
		if (!(w_N.w2 == 0 && w_N.w1 == 0 && w_N.w0 == 1))
			throw runtime_exception();
		w_N = w[N>>1];
		if (!(w_N.w2 == p.w2 && w_N.w1 == p.w1 && w_N.w0 == p.w0 - 1))
			throw runtime_exception();
	}
	#endif
	var bool squaring = (x == y);
	// Do an FFT of length N on x.
	{
		var sintL l;
		/* l = n-1 */ {
			var const uintC tmax = N>>1; // tmax = 2^(n-1)
			for (var uintC t = 0; t < tmax; t++) {
				var uintC i1 = t;
				var uintC i2 = i1 + tmax;
				// Butterfly: replace (x(i1),x(i2)) by
				// (x(i1) + x(i2), x(i1) - x(i2)).
				var fftp_word tmp;
				tmp = x[i2];
				subp(x[i1],tmp, x[i2]);
				addp(x[i1],tmp, x[i1]);
			}
		}
		for (l = n-2; l>=0; l--) {
			var const uintC smax = (uintC)1 << (n-1-l);
			var const uintC tmax = (uintC)1 << l;
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(n-1-l,s) << l;
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Butterfly: replace (x(i1),x(i2)) by
					// (x(i1) + w^exp*x(i2), x(i1) - w^exp*x(i2)).
					var fftp_word tmp;
					mulp(x[i2],w[exp], tmp);
					subp(x[i1],tmp, x[i2]);
					addp(x[i1],tmp, x[i1]);
				}
			}
		}
	}
	// Do an FFT of length N on y.
	if (!squaring) {
		var sintL l;
		/* l = n-1 */ {
			var uintC const tmax = N>>1; // tmax = 2^(n-1)
			for (var uintC t = 0; t < tmax; t++) {
				var uintC i1 = t;
				var uintC i2 = i1 + tmax;
				// Butterfly: replace (y(i1),y(i2)) by
				// (y(i1) + y(i2), y(i1) - y(i2)).
				var fftp_word tmp;
				tmp = y[i2];
				subp(y[i1],tmp, y[i2]);
				addp(y[i1],tmp, y[i1]);
			}
		}
		for (l = n-2; l>=0; l--) {
			var const uintC smax = (uintC)1 << (n-1-l);
			var const uintC tmax = (uintC)1 << l;
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(n-1-l,s) << l;
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Butterfly: replace (y(i1),y(i2)) by
					// (y(i1) + w^exp*y(i2), y(i1) - w^exp*y(i2)).
					var fftp_word tmp;
					mulp(y[i2],w[exp], tmp);
					subp(y[i1],tmp, y[i2]);
					addp(y[i1],tmp, y[i1]);
				}
			}
		}
	}
	// Multiply the transformed vectors into z.
	for (i = 0; i < N; i++)
		mulp(x[i],y[i], z[i]);
	// Undo an FFT of length N on z.
	{
		var uintL l;
		for (l = 0; l < n-1; l++) {
			var const uintC smax = (uintC)1 << (n-1-l);
			var const uintC tmax = (uintC)1 << l;
			#if FFTP_BACKWARD != CLEVER
			for (var uintC s = 0; s < smax; s++) {
				var uintC exp = bit_reverse(n-1-l,s) << l;
				#if FFTP_BACKWARD == RECIPROOT
				if (exp > 0)
					exp = N - exp; // negate exp (use w^-1 instead of w)
				#endif
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (z(i1),z(i2)) by
					// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/(2*w^exp)).
					var fftp_word sum;
					var fftp_word diff;
					addp(z[i1],z[i2], sum);
					subp(z[i1],z[i2], diff);
					shiftp(sum, z[i1]);
					mulp(diff,w[exp], diff); shiftp(diff, z[i2]);
				}
			}
			#else // FFTP_BACKWARD == CLEVER: clever handling of negative exponents
			/* s = 0, exp = 0 */ {
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (z(i1),z(i2)) by
					// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/(2*w^exp)),
					// with exp <-- 0.
					var fftp_word sum;
					var fftp_word diff;
					addp(z[i1],z[i2], sum);
					subp(z[i1],z[i2], diff);
					shiftp(sum, z[i1]);
					shiftp(diff, z[i2]);
				}
			}
			for (var uintC s = 1; s < smax; s++) {
				var uintC exp = bit_reverse(n-1-l,s) << l;
				exp = (N>>1) - exp; // negate exp (use w^-1 instead of w)
				for (var uintC t = 0; t < tmax; t++) {
					var uintC i1 = (s << (l+1)) + t;
					var uintC i2 = i1 + tmax;
					// Inverse Butterfly: replace (z(i1),z(i2)) by
					// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/(2*w^exp)),
					// with exp <-- (N/2 - exp).
					var fftp_word sum;
					var fftp_word diff;
					addp(z[i1],z[i2], sum);
					subp(z[i2],z[i1], diff); // note that w^(N/2) = -1
					shiftp(sum, z[i1]);
					mulp(diff,w[exp], diff); shiftp(diff, z[i2]);
				}
			}
			#endif
		}
		/* l = n-1 */ {
			var const uintC tmax = N>>1; // tmax = 2^(n-1)
			for (var uintC t = 0; t < tmax; t++) {
				var uintC i1 = t;
				var uintC i2 = i1 + tmax;
				// Inverse Butterfly: replace (z(i1),z(i2)) by
				// ((z(i1)+z(i2))/2, (z(i1)-z(i2))/2).
				var fftp_word sum;
				var fftp_word diff;
				addp(z[i1],z[i2], sum);
				subp(z[i1],z[i2], diff);
				shiftp(sum, z[i1]);
				shiftp(diff, z[i2]);
			}
		}
	}
	#if FFTP_BACKWARD == FORWARD
	// Swap z[i] and z[N-i] for 0 < i < N/2.
	for (i = (N>>1)-1; i > 0; i--) {
		var fftp_word tmp = z[i];
		z[i] = z[N-i];
		z[N-i] = tmp;
	}
	#endif
}

static void mulu_fft_modp (const uintD* sourceptr1, uintC len1,
                           const uintD* sourceptr2, uintC len2,
                           uintD* destptr)
// Es ist 2 <= len1 <= len2.
{
	// Methode:
	// source1 ist ein Stück der Länge N1, source2 ein oder mehrere Stücke
	// der Länge N2, mit N1+N2 <= N, wobei N Zweierpotenz ist.
	// sum(i=0..N-1, x_i b^i) * sum(i=0..N-1, y_i b^i) wird errechnet,
	// indem man die beiden Polynome
	// sum(i=0..N-1, x_i T^i), sum(i=0..N-1, y_i T^i)
	// multipliziert, und zwar durch Fourier-Transformation (s.o.).
	var uint32 n;
	integerlengthC(len1-1, n=); // 2^(n-1) < len1 <= 2^n
	var uintC len = (uintC)1 << n; // kleinste Zweierpotenz >= len1
	// Wählt man N = len, so hat man ceiling(len2/(len-len1+1)) * FFT(len).
	// Wählt man N = 2*len, so hat man ceiling(len2/(2*len-len1+1)) * FFT(2*len).
	// Wir wählen das billigere von beiden:
	// Bei ceiling(len2/(len-len1+1)) <= 2 * ceiling(len2/(2*len-len1+1))
	// nimmt man N = len, bei ....... > ........ dagegen N = 2*len.
	// (Wahl von N = 4*len oder mehr bringt nur in Extremfällen etwas.)
	if (len2 > 2 * (len-len1+1) * (len2 <= (2*len-len1+1) ? 1 : ceiling(len2,(2*len-len1+1)))) {
		n = n+1;
		len = len << 1;
	}
	var const uintC N = len; // N = 2^n
	CL_ALLOCA_STACK;
	var fftp_word* const x = cl_alloc_array(fftp_word,N);
	var fftp_word* const y = cl_alloc_array(fftp_word,N);
	#ifdef DEBUG_FFTP
	var fftp_word* const z = cl_alloc_array(fftp_word,N);
	#else
	var fftp_word* const z = x; // put z in place of x - saves memory
	#endif
	var uintD* const tmpprod = cl_alloc_array(uintD,len1+1);
	var uintP i;
	var uintC destlen = len1+len2;
	clear_loop_lsp(destptr,destlen);
	do {
		var uintC len2p; // length of a piece of source2
		len2p = N - len1 + 1;
		if (len2p > len2)
			len2p = len2;
		// len2p = min(N-len1+1,len2).
		if (len2p == 1) {
			// cheap case
			var uintD* tmpptr = arrayLSDptr(tmpprod,len1+1);
			mulu_loop_lsp(lspref(sourceptr2,0),sourceptr1,tmpptr,len1);
			if (addto_loop_lsp(tmpptr,destptr,len1+1))
				if (inc_loop_lsp(destptr lspop (len1+1),destlen-(len1+1)))
					throw runtime_exception();
		} else {
			var uintC destlenp = len1 + len2p - 1;
			// destlenp = min(N,destlen-1).
			var bool squaring = ((sourceptr1 == sourceptr2) && (len1 == len2p));
			// Fill factor x.
			{
				for (i = 0; i < len1; i++) {
					x[i].w0 = lspref(sourceptr1,i);
					x[i].w1 = 0;
					x[i].w2 = 0;
				}
				for (i = len1; i < N; i++) {
					x[i].w0 = 0;
					x[i].w1 = 0;
					x[i].w2 = 0;
				}
			}
			// Fill factor y.
			if (!squaring) {
				for (i = 0; i < len2p; i++) {
					y[i].w0 = lspref(sourceptr2,i);
					y[i].w1 = 0;
					y[i].w2 = 0;
				}
				for (i = len2p; i < N; i++) {
					y[i].w0 = 0;
					y[i].w1 = 0;
					y[i].w2 = 0;
				}
			}
			// Multiply.
			if (!squaring)
				fftp_convolution(n,N, &x[0], &y[0], &z[0]);
			else
				fftp_convolution(n,N, &x[0], &x[0], &z[0]);
			#ifdef DEBUG_FFTP
			// Check result.
			for (i = 0; i < N; i++)
				if (!(z[i].w2 < N))
					throw runtime_exception();
			#endif
			// Add result to destptr[-destlen..-1]:
			{
				var uintD* ptr = destptr;
				// ac2|ac1|ac0 are an accumulator.
				var uint32 ac0 = 0;
				var uint32 ac1 = 0;
				var uint32 ac2 = 0;
				var uint32 tmp;
				for (i = 0; i < destlenp; i++) {
					// Add z[i] to the accumulator.
					tmp = z[i].w0;
					if ((ac0 += tmp) < tmp) {
						if (++ac1 == 0)
							++ac2;
					}
					tmp = z[i].w1;
					if ((ac1 += tmp) < tmp)
						++ac2;
					tmp = z[i].w2;
					ac2 += tmp;
					// Add the accumulator's least significant word to destptr:
					tmp = lspref(ptr,0);
					if ((ac0 += tmp) < tmp) {
						if (++ac1 == 0)
							++ac2;
					}
					lspref(ptr,0) = ac0;
					lsshrink(ptr);
					ac0 = ac1;
					ac1 = ac2;
					ac2 = 0;
				}
				// ac2 = 0.
				if (ac1 > 0) {
					if (!((i += 2) <= destlen))
						throw runtime_exception();
					tmp = lspref(ptr,0);
					if ((ac0 += tmp) < tmp)
						++ac1;
					lspref(ptr,0) = ac0;
					lsshrink(ptr);
					tmp = lspref(ptr,0);
					ac1 += tmp;
					lspref(ptr,0) = ac1;
					lsshrink(ptr);
					if (ac1 < tmp)
						if (inc_loop_lsp(ptr,destlen-i))
							throw runtime_exception();
				} else if (ac0 > 0) {
					if (!((i += 1) <= destlen))
						throw runtime_exception();
					tmp = lspref(ptr,0);
					ac0 += tmp;
					lspref(ptr,0) = ac0;
					lsshrink(ptr);
					if (ac0 < tmp)
						if (inc_loop_lsp(ptr,destlen-i))
							throw runtime_exception();
				}
			}
			#ifdef DEBUG_FFTP
			// If destlenp < N, check that the remaining z[i] are 0.
			for (i = destlenp; i < N; i++)
				if (z[i].w2 > 0 || z[i].w1 > 0 || z[i].w0 > 0)
					throw runtime_exception();
			#endif
		}
		// Decrement len2.
		destptr = destptr lspop len2p;
		destlen -= len2p;
		sourceptr2 = sourceptr2 lspop len2p;
		len2 -= len2p;
	} while (len2 > 0);
}

#undef FFT_P_94
#undef FFT_P_92
#undef w0
#undef w1
#undef w2
#undef W3
