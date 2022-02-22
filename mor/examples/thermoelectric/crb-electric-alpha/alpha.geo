h=0.3;

ZB=0.;
ZT=5.;
Zs=0.5;
Ze=4.5;
Zx=0.05;
RI=0.5;
RE=1.;
nbTour=2;
n=20;

ZBox = 2.5;
RBox = 0.2;

Macro cyl
    cylP[] = {};
    cylL[] = {};
    cylS[] = {};
    cylV[] = {};
    For i In { 0 : n }
    	angle = i*2*Pi/n;

	If ( i != n )
	    cylP[] += {newp}; Point(cylP[4*i]) = {Ri*Cos(angle),Ri*Sin(angle),Z,h};
    	    cylP[] += {newp}; Point(cylP[4*i+1]) = {Re*Cos(angle),Re*Sin(angle),Z,h};
	    cylP[] += {newp}; Point(cylP[4*i+2]) = {Ri*Cos(angle),Ri*Sin(angle),Z+L,h};
    	    cylP[] += {newp}; Point(cylP[4*i+3]) = {Re*Cos(angle),Re*Sin(angle),Z+L,h};

	    cylL[] += {newl}; Line(cylL[8*i+0]) = {cylP[4*i+0],cylP[4*i+1]};
            cylL[] += {newl}; Line(cylL[8*i+1]) = {cylP[4*i+2],cylP[4*i+3]};
            cylL[] += {newl}; Line(cylL[8*i+2]) = {cylP[4*i+0],cylP[4*i+2]};
            cylL[] += {newl}; Line(cylL[8*i+3]) = {cylP[4*i+1],cylP[4*i+3]};
	Else
	    cylL[] += {newl};
	    cylL[] += {newl};
	    cylL[] += {newl};
	    cylL[] += {newl};
	EndIf
	
	If ( i == 0 )
	    cylL[] += {newl};
	    cylL[] += {newl};
	    cylL[] += {newl};
	    cylL[] += {newl};
    	Else
    	    cylL[] += {newl}; Line(cylL[8*i+4]) = {cylP[4*(i-1)+0],cylP[(4*i+0)%(4*n)]};
    	    cylL[] += {newl}; Line(cylL[8*i+5]) = {cylP[4*(i-1)+1],cylP[(4*i+1)%(4*n)]};
    	    cylL[] += {newl}; Line(cylL[8*i+6]) = {cylP[4*(i-1)+2],cylP[(4*i+2)%(4*n)]};
    	    cylL[] += {newl}; Line(cylL[8*i+7]) = {cylP[4*(i-1)+3],cylP[(4*i+3)%(4*n)]};
	EndIf


	If ( i == 0 )
	    cylS[] += {news};
	    cylS[] += {news};
	    cylS[] += {news};
	    cylS[] += {news};
	Else
	    ll = newll; Line Loop(ll) = {cylL[(8*i)%(8*n)],-cylL[8*i+5],-cylL[8*(i-1)],cylL[8*i+4]};
	    cylS[] += {news}; Plane Surface(cylS[5*i]) = {ll};
	    ll = newll; Line Loop(ll) = {cylL[(8*i+1)%(8*n)],-cylL[8*i+7],-cylL[8*(i-1)+1],cylL[8*i+6]};
	    cylS[] += {news}; Plane Surface(cylS[5*i+1]) = {ll};
	    ll = newll; Line Loop(ll) = {cylL[(8*i+2)%(8*n)],-cylL[8*i+6],-cylL[8*(i-1)+2],cylL[8*i+4]};
	    cylS[] += {news}; Plane Surface(cylS[5*i+2]) = {ll};
	    ll = newll; Line Loop(ll) = {cylL[(8*i+3)%(8*n)],-cylL[8*i+7],-cylL[8*(i-1)+3],cylL[8*i+5]};
	    cylS[] += {news}; Plane Surface(cylS[5*i+3]) = {ll};
	EndIf

	If ( i != n )
	    ll = newll; Line Loop(ll) = {cylL[8*i+0],cylL[8*i+3],-cylL[8*i+1],-cylL[8*i+2]};
	    cylS[] += {news}; Plane Surface(cylS[5*i+4]) = {ll};
	EndIf

	If ( i != 0 )
	    sl = newsl; Surface Loop(sl) = {cylS[5*i+0],cylS[5*i+1],cylS[5*i+2],cylS[5*i+3],cylS[(5*i+4)%(5*n)],cylS[5*(i-1)+4]};
	    cylV[] += {newv}; Volume(cylV[i-1]) = {sl};
	EndIf
    EndFor
Return

Macro Helix
    helP[] = {};
    helL[] = {};
    helS[] = {};
    helV[] = {};
    N = nbTour*n;
    For i In { 0 : N }
    	angle = i*2*Pi/n;

	// Points
	If ( i == 0 )
   	    helP[] += {baseIP};
	    helP[] += {baseEP};
    	    helP[] += {newp}; Point(helP[4*i+2]) = {Ri*Cos(angle),Ri*Sin(angle),Zs+i*(Ze-Zs-Zx)/N+Zx,h};
    	    helP[] += {newp}; Point(helP[4*i+3]) = {Re*Cos(angle),Re*Sin(angle),Zs+i*(Ze-Zs-Zx)/N+Zx,h};
	ElseIf ( i == N )
	    helP[] += {newp}; Point(helP[4*i]) = {Ri*Cos(angle),Ri*Sin(angle),Zs+i*(Ze-Zs-Zx)/N,h};
    	    helP[] += {newp}; Point(helP[4*i+1]) = {Re*Cos(angle),Re*Sin(angle),Zs+i*(Ze-Zs-Zx)/N,h};
   	    helP[] += {topIP};
	    helP[] += {topEP};
	Else
	    helP[] += {newp}; Point(helP[4*i]) = {Ri*Cos(angle),Ri*Sin(angle),Zs+i*(Ze-Zs-Zx)/N,h};
    	    helP[] += {newp}; Point(helP[4*i+1]) = {Re*Cos(angle),Re*Sin(angle),Zs+i*(Ze-Zs-Zx)/N,h};
    	    helP[] += {newp}; Point(helP[4*i+2]) = {Ri*Cos(angle),Ri*Sin(angle),Zs+i*(Ze-Zs-Zx)/N+Zx,h};
    	    helP[] += {newp}; Point(helP[4*i+3]) = {Re*Cos(angle),Re*Sin(angle),Zs+i*(Ze-Zs-Zx)/N+Zx,h};
	EndIf

    	// Lines transverse
	If( i == 0 )
	    helL[] += {baseLine};
    	    helL[] += {newl}; Line(helL[8*i+1]) = {helP[4*i+2], helP[4*i+3]};
	ElseIf ( i == N )
    	    helL[] += {newl}; Line(helL[8*i]) = {helP[4*i],helP[4*i+1]};
	    helL[] += {topLine};
	Else
	    helL[] += {newl}; Line(helL[8*i]) = {helP[4*i],helP[4*i+1]};
	    helL[] += {newl}; Line(helL[8*i+1]) = {helP[4*i+2], helP[4*i+3]};
	EndIf

    	// Lines between helix
    	helL[] += {newl}; Line(helL[8*i+2]) = {helP[4*i],helP[4*i+2]};
    	helL[] += {newl}; Line(helL[8*i+3]) = {helP[4*i+1],helP[4*i+3]};

	// Lines across the helix
    	If ( i != 0 )
            helL[] += {newl}; Line(helL[8*i+4]) = {helP[4*i],helP[4*(i-1)]};
	    helL[] += {newl}; Line(helL[8*i+5]) = {helP[4*i+1],helP[4*(i-1)+1]};
	    helL[] += {newl}; Line(helL[8*i+6]) = {helP[4*i+2],helP[4*(i-1)+2]};
	    helL[] += {newl}; Line(helL[8*i+7]) = {helP[4*i+3],helP[4*(i-1)+3]};
	Else
	    helL[] += {newl};
	    helL[] += {newl};
	    helL[] += {newl};
	    helL[] += {newl};
    	EndIf

        ll = newll; Line Loop(ll) = {helL[8*i+2],helL[8*i+1],-helL[8*i+3],-helL[8*i]};
        helS[] += {news}; Ruled Surface(helS[5*i]) = {ll};

        If ( i != 0 )
            ll = newll; Line Loop(ll) = {helL[8*i],helL[8*i+5],-helL[8*(i-1)],-helL[8*i+4]};
            helS[] += {newll}; Ruled Surface(helS[5*i+1]) = {ll};
            ll = newll; Line Loop(ll) = {helL[8*i+1],helL[8*i+7],-helL[8*(i-1)+1],-helL[8*i+6]};
            helS[] += {newll}; Ruled Surface(helS[5*i+2]) = {ll};
    
	    ll = newll; Line Loop(ll) = {helL[8*i+2],helL[8*i+6],-helL[8*(i-1)+2],-helL[8*i+4]};
            helS[] += {news}; Ruled Surface(helS[5*i+3]) = {ll};
            ll = newll; Line Loop(ll) = {helL[8*i+3],helL[8*i+7],-helL[8*(i-1)+3],-helL[8*i+5]};
            helS[] += {news}; Ruled Surface(helS[5*i+4]) = {ll};
	Else
	    ll = newll; helS[] += {news};
	    ll = newll; helS[] += {news};
	    ll = newll; helS[] += {news};
	    ll = newll; helS[] += {news};
        EndIf

	If( i != 0 )
    	    sl = newsl; Surface Loop(sl) = {helS[5*(i-1)],helS[5*i],helS[5*i+1],helS[5*i+2],helS[5*i+3],helS[5*i+4]};
    	    helV[] += {newv}; Volume(helV[i-1]) = {sl};
    	EndIf
    EndFor
Return

Macro CondDown
    condDL[] = {};
    condDS[] = {};
    condDV[] = {};
    N = n*nbTour;
    For i In { 1 : n/2-1 }
    	condDL[] += {newl}; Line(condDL[4*(i-1)+0]) = {baseP[4*i+2],hel1P[4*i]};
    	condDL[] += {newl}; Line(condDL[4*(i-1)+1]) = {baseP[4*i+3],hel1P[4*i+1]};
    	condDL[] += {newl}; Line(condDL[4*(i-1)+2]) = {baseP[4*(n/2+i)+2],hel2P[4*i]};
    	condDL[] += {newl}; Line(condDL[4*(i-1)+3]) = {baseP[4*(n/2+i)+3],hel2P[4*i+1]};

	ll = newll; Line Loop(ll) = {baseL[8*i+1],condDL[4*(i-1)+1],-hel1L[8*i],-condDL[4*(i-1)]};
	condDS[] += {news}; Ruled Surface(condDS[6*(i-1)]) = {ll};
	ll = newll; Line Loop(ll) = {baseL[8*(n/2+i)+1],condDL[4*(i-1)+3],-hel2L[8*i],-condDL[4*(i-1)+2]};
	condDS[] += {news}; Ruled Surface(condDS[6*(i-1)+1]) = {ll};

	lll1[] = {baseL[8*i+6],condDL[4*(i-1)],hel1L[8*i+4]};
	lll2[] = {baseL[8*(n/2+i)+6],condDL[4*(i-1)+2],hel2L[8*i+4]};
	lll3[] = {baseL[8*i+7],condDL[4*(i-1)+1],hel1L[8*i+5]};
	lll4[] = {baseL[8*(n/2+i)+7],condDL[4*(i-1)+3],hel2L[8*i+5]};
	If ( i != 1 )
	   lll1[] += {-condDL[4*(i-2)]};
	   lll2[] += {-condDL[4*(i-2)+2]};
	   lll3[] += {-condDL[4*(i-2)+1]};
	   lll4[] += {-condDL[4*(i-2)+3]};
	EndIf
	ll = newll; Line Loop(ll) = {lll1[]};
	condDS[] += {news}; Ruled Surface(condDS[6*(i-1)+2]) = {ll};
	ll = newll; Line Loop(ll) = {lll2[]};
	condDS[] += {news}; Ruled Surface(condDS[6*(i-1)+3]) = {ll};
	ll = newll; Line Loop(ll) = {lll3[]};
	condDS[] += {news}; Ruled Surface(condDS[6*(i-1)+4]) = {ll};
	ll = newll; Line Loop(ll) = {lll4[]};
	condDS[] += {news}; Ruled Surface(condDS[6*(i-1)+5]) = {ll};

	ssl1[] = {condDS[6*(i-1)],condDS[6*(i-1)+2],condDS[6*(i-1)+4],baseS[5*i+1],hel1S[5*i+1]};
	ssl2[] = {condDS[6*(i-1)+1],condDS[6*(i-1)+3],condDS[6*(i-1)+5],baseS[5*(n/2+i)+1],hel2S[5*i+1]};
	If ( i != 1 )
	   ssl1[] += {condDS[6*(i-2)]};
	   ssl2[] += {condDS[6*(i-2)+1]};
	EndIf
	sl = newsl; Surface Loop(sl) = {ssl1[]};
	condDV[] += {newv}; Volume(condDV[2*(i-1)]) = {sl};
	sl = newsl; Surface Loop(sl) = {ssl2[]};
	condDV[] += {newv}; Volume(condDV[2*(i-1)+1]) = {sl};
    EndFor

    condDL[] += {newl}; Line(condDL[4*(n/2-1)+0]) = {hel1P[2],baseP[4*(n-1)+2]};
    condDL[] += {newl}; Line(condDL[4*(n/2-1)+1]) = {hel1P[3],baseP[4*(n-1)+3]};
    condDL[] += {newl}; Line(condDL[4*(n/2-1)+2]) = {hel2P[2],baseP[4*(n/2-1)+2]};
    condDL[] += {newl}; Line(condDL[4*(n/2-1)+3]) = {hel2P[3],baseP[4*(n/2-1)+3]};
    condDL[] += {newl}; Line(condDL[4*(n/2)+0]) = {hel1P[2],hel2P[4*(n/2)]};
    condDL[] += {newl}; Line(condDL[4*(n/2)+1]) = {hel1P[3],hel2P[4*(n/2)+1]};
    condDL[] += {newl}; Line(condDL[4*(n/2)+2]) = {hel2P[2],hel1P[4*(n/2)]};
    condDL[] += {newl}; Line(condDL[4*(n/2)+3]) = {hel2P[3],hel1P[4*(n/2)+1]};

    ll = newll; Line Loop(ll) = {hel1L[1],condDL[4*(n/2-1)+1],-baseL[8*(n-1)+1],-condDL[4*(n/2-1)+0]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2-1)]) = {ll};
    ll = newll; Line Loop(ll) = {hel2L[1],condDL[4*(n/2-1)+3],-baseL[8*(n/2-1)+1],-condDL[4*(n/2-1)+2]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2-1)+1]) = {ll};
    ll = newll; Line Loop(ll) = {hel1L[2],condDL[4*(n/2-1)],baseL[8*n+6]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2-1)+2]) = {ll};
    ll = newll; Line Loop(ll) = {hel2L[2],condDL[4*(n/2-1)+2],baseL[8*n/2+6]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2-1)+3]) = {ll};
    ll = newll; Line Loop(ll) = {hel1L[3],condDL[4*(n/2-1)+1],baseL[8*n+7]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2-1)+4]) = {ll};
    ll = newll; Line Loop(ll) = {hel2L[3],condDL[4*(n/2-1)+3],baseL[8*n/2+7]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2-1)+5]) = {ll};

    sl = newsl; Surface Loop(sl) = {condDS[6*(n/2-1)],condDS[6*(n/2-1)+2],condDS[6*(n/2-1)+4],baseS[5*n+1],hel1S[0]};
    condDV[] += {newv}; Volume(condDV[2*(n/2-1)]) = {sl};
    sl = newsl; Surface Loop(sl) = {condDS[6*(n/2-1)+1],condDS[6*(n/2-1)+3],condDS[6*(n/2-1)+5],baseS[5*(n/2)+1],hel2S[0]};
    condDV[] += {newv}; Volume(condDV[2*(n/2-1)+1]) = {sl};

    ll = newll; Line Loop(ll) = {hel1L[1],condDL[4*(n/2)+1],-hel2L[8*(n/2)],-condDL[4*(n/2)]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2)]) = {ll};
    ll = newll; Line Loop(ll) = {hel2L[1],condDL[4*(n/2)+3],-hel1L[8*(n/2)],-condDL[4*(n/2)+2]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2)+1]) = {ll};
    ll = newll; Line Loop(ll) = {condDL[4*(n/2)],hel2L[8*(n/2)+4],-condDL[4*(n/2-2)+2],-condDL[4*(n/2-1)]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2)+2]) = {ll};
    ll = newll; Line Loop(ll) = {condDL[4*(n/2)+2],hel1L[8*(n/2)+4],-condDL[4*(n/2-2)],-condDL[4*(n/2-1)+2]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2)+3]) = {ll};
    ll = newll; Line Loop(ll) = {condDL[4*(n/2)+1],hel2L[8*(n/2)+5],-condDL[4*(n/2-2)+3],-condDL[4*(n/2-1)+1]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2)+4]) = {ll};
    ll = newll; Line Loop(ll) = {condDL[4*(n/2)+3],hel1L[8*(n/2)+5],-condDL[4*(n/2-2)+1],-condDL[4*(n/2-1)+3]};
    condDS[] += {news}; Ruled Surface(condDS[6*(n/2)+5]) = {ll};

    sl = newsl; Surface Loop(sl) = {condDS[6*(n/2)],condDS[6*(n/2)+2],condDS[6*(n/2)+4],condDS[6*(n/2-1)],hel2S[5*(n/2)+1],condDS[6*(n/2-2)+1]};
    condDV[] += {newv}; Volume(condDV[2*(n/2)]) = {sl};
    sl = newsl; Surface Loop(sl) = {condDS[6*(n/2)+1],condDS[6*(n/2)+3],condDS[6*(n/2)+5],condDS[6*(n/2-1)+1],hel1S[5*(n/2)+1],condDS[6*(n/2-2)]};
    condDV[] += {newv}; Volume(condDV[2*(n/2)+1]) = {sl};    
Return

Macro CondUp
    condUL[] = {};
    condUS[] = {};
    condUV[] = {};
    N = n*nbTour;
    For i In { 1 : n/2-1 }
    	condUL[] += {newl}; Line(condUL[4*(i-1)+0]) = {hel1P[4*(N-i)+2],topP[4*(n-i)+0]};
    	condUL[] += {newl}; Line(condUL[4*(i-1)+1]) = {hel1P[4*(N-i)+3],topP[4*(n-i)+1]};
    	condUL[] += {newl}; Line(condUL[4*(i-1)+2]) = {hel2P[4*(N-i)+2],topP[4*(n/2-i)+0]};
    	condUL[] += {newl}; Line(condUL[4*(i-1)+3]) = {hel2P[4*(N-i)+3],topP[4*(n/2-i)+1]};
    	
	ll = newll; Line Loop(ll) = {hel1L[8*(N-i)+1],condUL[4*(i-1)+1],-topL[8*(n-i)],-condUL[4*(i-1)]};
	condUS[] += {news}; Ruled Surface(condUS[6*(i-1)]) = {ll};
	ll = newll; Line Loop(ll) = {hel2L[8*(N-i)+1],condUL[4*(i-1)+3],-topL[8*(n/2-i)],-condUL[4*(i-1)+2]};
	condUS[] += {news}; Ruled Surface(condUS[6*(i-1)+1]) = {ll};

	lll1[] = {-topL[8*(n-i+1)+4],-condUL[4*(i-1)+0],-hel1L[8*(N-i+1)+6]};
	lll2[] = {-topL[8*(n/2-i+1)+4],-condUL[4*(i-1)+2],-hel2L[8*(N-i+1)+6]};
	lll3[] = {-topL[8*(n-i+1)+5],-condUL[4*(i-1)+1],-hel1L[8*(N-i+1)+7]};
	lll4[] = {-topL[8*(n/2-i+1)+5],-condUL[4*(i-1)+3],-hel2L[8*(N-i+1)+7]};
	If ( i != 1 )
	   lll1[] += {condUL[4*(i-2)]};
	   lll2[] += {condUL[4*(i-2)+2]};
	   lll3[] += {condUL[4*(i-2)+1]};
	   lll4[] += {condUL[4*(i-2)+3]};
	EndIf
	ll = newll; Line Loop(ll) = {lll1[]};
	condUS[] += {news}; Ruled Surface(condUS[6*(i-1)+2]) = {ll};
	ll = newll; Line Loop(ll) = {lll2[]};
	condUS[] += {news}; Ruled Surface(condUS[6*(i-1)+3]) = {ll};
	ll = newll; Line Loop(ll) = {lll3[]};
	condUS[] += {news}; Ruled Surface(condUS[6*(i-1)+4]) = {ll};
	ll = newll; Line Loop(ll) = {lll4[]};
	condUS[] += {news}; Ruled Surface(condUS[6*(i-1)+5]) = {ll};

	ssl1[] = {condUS[6*(i-1)],condUS[6*(i-1)+2],condUS[6*(i-1)+4],topS[5*(n-i+1)],hel1S[5*(N-i+1)+2]};
	ssl2[] = {condUS[6*(i-1)+1],condUS[6*(i-1)+3],condUS[6*(i-1)+5],topS[5*(n/2-i+1)],hel2S[5*(N-i+1)+2]};
	If ( i != 1 )
	   ssl1[] += {condUS[6*(i-2)]};
	   ssl2[] += {condUS[6*(i-2)+1]};
	EndIf
	sl = newsl; Surface Loop(sl) = {ssl1[]};
	condUV[] += {newv}; Volume(condUV[2*(i-1)]) = {sl};
	sl = newsl; Surface Loop(sl) = {ssl2[]};
	condUV[] += {newv}; Volume(condUV[2*(i-1)+1]) = {sl};
    EndFor

    condUL[] += {newl}; Line(condUL[4*(n/2-1)+0]) = {hel1P[4*N],topP[4]};
    condUL[] += {newl}; Line(condUL[4*(n/2-1)+1]) = {hel1P[4*N+1],topP[5]};
    condUL[] += {newl}; Line(condUL[4*(n/2-1)+2]) = {hel2P[4*N],topP[4*(n/2+1)]};
    condUL[] += {newl}; Line(condUL[4*(n/2-1)+3]) = {hel2P[4*N+1],topP[4*(n/2+1)+1]};
    condUL[] += {newl}; Line(condUL[4*(n/2)+0]) = {hel2P[4*(N-n/2)+2],hel1P[4*N]};
    condUL[] += {newl}; Line(condUL[4*(n/2)+1]) = {hel2P[4*(N-n/2)+3],hel1P[4*N+1]};
    condUL[] += {newl}; Line(condUL[4*(n/2)+2]) = {hel1P[4*(N-n/2)+2],hel2P[4*N]};
    condUL[] += {newl}; Line(condUL[4*(n/2)+3]) = {hel1P[4*(N-n/2)+3],hel2P[4*N+1]};

    ll = newll; Line Loop(ll) = {topL[8],-condUL[4*(n/2-1)+1],-hel1L[8*N],condUL[4*(n/2-1)+0]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2-1)]) = {ll};
    ll = newll; Line Loop(ll) = {topL[8*(n/2+1)],-condUL[4*(n/2-1)+3],-hel2L[8*N],condUL[4*(n/2-1)+2]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2-1)+1]) = {ll};
    ll = newll; Line Loop(ll) = {hel1L[8*N+2],topL[12],-condUL[4*(n/2-1)]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2-1)+2]) = {ll};
    ll = newll; Line Loop(ll) = {hel2L[8*N+2],topL[8*(n/2+1)+4],-condUL[4*(n/2-1)+2]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2-1)+3]) = {ll};
    ll = newll; Line Loop(ll) = {hel1L[8*N+3],topL[13],-condUL[4*(n/2-1)+1]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2-1)+4]) = {ll};
    ll = newll; Line Loop(ll) = {hel2L[8*N+3],topL[8*(n/2+1)+5],-condUL[4*(n/2-1)+3]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2-1)+5]) = {ll};

    sl = newsl; Surface Loop(sl) = {condUS[6*(n/2-1)],condUS[6*(n/2-1)+2],condUS[6*(n/2-1)+4],topS[5],hel1S[5*N]};
    condUV[] += {newv}; Volume(condUV[2*(n/2-1)]) = {sl};
    sl = newsl; Surface Loop(sl) = {condUS[6*(n/2-1)+1],condUS[6*(n/2-1)+3],condUS[6*(n/2-1)+5],topS[5*(n/2+1)],hel2S[5*N]};
    condUV[] += {newv}; Volume(condUV[2*(n/2-1)+1]) = {sl};

    ll = newll; Line Loop(ll) = {hel1L[8*N],-condUL[4*(n/2)+1],-hel2L[8*(N-n/2)+1],condUL[4*(n/2)]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2)]) = {ll};
    ll = newll; Line Loop(ll) = {hel2L[8*N],-condUL[4*(n/2)+3],-hel1L[8*(N-n/2)+1],condUL[4*(n/2)+2]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2)+1]) = {ll};
    ll = newll; Line Loop(ll) = {condUL[4*(n/2)],condUL[4*(n/2-1)],-condUL[4*(n/2-2)+2],hel2L[8*(N-n/2+1)+6]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2)+2]) = {ll};
    ll = newll; Line Loop(ll) = {condUL[4*(n/2)+2],condUL[4*(n/2-1)+2],-condUL[4*(n/2-2)],hel1L[8*(N-n/2+1)+6]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2)+3]) = {ll};
    ll = newll; Line Loop(ll) = {condUL[4*(n/2)+1],condUL[4*(n/2-1)+1],-condUL[4*(n/2-2)+3],hel2L[8*(N-n/2+1)+7]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2)+4]) = {ll};
    ll = newll; Line Loop(ll) = {condUL[4*(n/2)+3],condUL[4*(n/2-1)+3],-condUL[4*(n/2-2)+1],hel1L[8*(N-n/2+1)+7]};
    condUS[] += {news}; Ruled Surface(condUS[6*(n/2)+5]) = {ll};

    sl = newsl; Surface Loop(sl) = {condUS[6*(n/2)],condUS[6*(n/2)+2],condUS[6*(n/2)+4],condUS[6*(n/2-1)],hel2S[5*(N-n/2+1)+2],condUS[6*(n/2-2)+1]};
    condUV[] += {newv}; Volume(condUV[2*(n/2)]) = {sl};
    sl = newsl; Surface Loop(sl) = {condUS[6*(n/2)+1],condUS[6*(n/2)+3],condUS[6*(n/2)+5],condUS[6*(n/2-1)+1],hel1S[5*(N-n/2+1)+2],condUS[6*(n/2-2)]};
    condUV[] += {newv}; Volume(condUV[2*(n/2)+1]) = {sl};    
Return

Macro CondBet
    condBL[] = {};
    condBS[] = {};
    condBV[] = {};
    N = n*nbTour;
    condBL[] += {condDL[4*(n/2)+0],condDL[4*(n/2)+1],condDL[4*(n/2)+2],condDL[4*(n/2)+3]};
    condBS[] += {condDS[6*(n/2)],condDS[6*(n/2)+1],news,news,news,news};
    For i In { 1 : N-n/2 }
    	If ( i == N-n/2 )
	   condBL[] += {condUL[4*(n/2)+2]};
	   condBL[] += {condUL[4*(n/2)+3]};
	   condBL[] += {condUL[4*(n/2)+0]};
	   condBL[] += {condUL[4*(n/2)+1]};
	   condBS[] += {condUS[6*(n/2)+1]};
	   condBS[] += {condUS[6*(n/2)+0]};
	Else
	   condBL[] += {newl}; Line(condBL[4*i]) = {hel1P[4*i+2],hel2P[4*(n/2+i)]};
	   condBL[] += {newl}; Line(condBL[4*i+1]) = {hel1P[4*i+3],hel2P[4*(n/2+i)+1]};
	   condBL[] += {newl}; Line(condBL[4*i+2]) = {hel2P[4*i+2],hel1P[4*(n/2+i)]};
	   condBL[] += {newl}; Line(condBL[4*i+3]) = {hel2P[4*i+3],hel1P[4*(n/2+i)+1]};

	   ll = newll; Line Loop(ll) = {condBL[4*i],hel2L[8*(n/2+i)],-condBL[4*i+1],-hel1L[8*i+1]};
	   condBS[] += {news}; Ruled Surface(condBS[6*i]) = {ll};
	   ll = newll; Line Loop(ll) = {condBL[4*i+2],hel1L[8*(n/2+i)],-condBL[4*i+3],-hel2L[8*i+1]};
	   condBS[] += {news}; Ruled Surface(condBS[6*i+1]) = {ll};
	EndIf

	ll = newll; Line Loop(ll) = {condBL[4*(i-1)],-hel2L[8*(n/2+i)+4],-condBL[4*i],hel1L[8*i+6]};
	condBS[] += {news}; Ruled Surface(condBS[6*i+2]) = {ll};
	ll = newll; Line Loop(ll) = {condBL[4*(i-1)+2],-hel1L[8*(n/2+i)+4],-condBL[4*i+2],hel2L[8*i+6]};
	condBS[] += {news}; Ruled Surface(condBS[6*i+3]) = {ll};

	ll = newll; Line Loop(ll) = {condBL[4*(i-1)+1],-hel2L[8*(n/2+i)+5],-condBL[4*i+1],hel1L[8*i+7]};
	condBS[] += {news}; Ruled Surface(condBS[6*i+4]) = {ll};
	ll = newll; Line Loop(ll) = {condBL[4*(i-1)+3],-hel1L[8*(n/2+i)+5],-condBL[4*i+3],hel2L[8*i+7]};
	condBS[] += {news}; Ruled Surface(condBS[6*i+5]) = {ll};

	sl = newsl; Surface Loop(sl) = {condBS[6*i],condBS[6*i+2],condBS[6*i+4],condBS[6*(i-1)],hel1S[5*i+2],hel2S[5*(n/2+i)+1]};
	condBV[] += {newv}; Volume(condBV[2*(i-1)]) = {sl};

	sl = newsl; Surface Loop(sl) = {condBS[6*i+1],condBS[6*i+3],condBS[6*i+5],condBS[6*(i-1)+1],hel2S[5*i+2],hel1S[5*(n/2+i)+1]};
	condBV[] += {newv}; Volume(condBV[2*(i-1)+1]) = {sl};
    EndFor
Return

Macro Box
    boxP[] = {};
    boxL[] = {};
    boxS[] = {};
    boxV = 0;

    boxP[] += {newp}; Point(boxP[0]) = {0,0,ZBox,h};
    boxP[] += {newp}; Point(boxP[1]) = {RBox,0,ZBox,h};
    boxP[] += {newp}; Point(boxP[2]) = {0,RBox,ZBox,h};
    boxP[] += {newp}; Point(boxP[3]) = {-RBox,0,ZBox,h};
    boxP[] += {newp}; Point(boxP[4]) = {0,-RBox,ZBox,h};
    boxP[] += {newp}; Point(boxP[5]) = {0,0,ZBox+RBox,h};
    boxP[] += {newp}; Point(boxP[6]) = {0,0,ZBox-RBox,h};

    boxL[] += {newl}; Circle(boxL[0]) = {boxP[1],boxP[0],boxP[2]};
    boxL[] += {newl}; Circle(boxL[1]) = {boxP[2],boxP[0],boxP[3]};
    boxL[] += {newl}; Circle(boxL[2]) = {boxP[3],boxP[0],boxP[4]};
    boxL[] += {newl}; Circle(boxL[3]) = {boxP[4],boxP[0],boxP[1]};
    boxL[] += {newl}; Circle(boxL[4]) = {boxP[1],boxP[0],boxP[5]};
    boxL[] += {newl}; Circle(boxL[5]) = {boxP[1],boxP[0],boxP[6]};
    boxL[] += {newl}; Circle(boxL[6]) = {boxP[2],boxP[0],boxP[5]};
    boxL[] += {newl}; Circle(boxL[7]) = {boxP[2],boxP[0],boxP[6]};
    boxL[] += {newl}; Circle(boxL[8]) = {boxP[3],boxP[0],boxP[5]};
    boxL[] += {newl}; Circle(boxL[9]) = {boxP[3],boxP[0],boxP[6]};
    boxL[] += {newl}; Circle(boxL[10]) = {boxP[4],boxP[0],boxP[5]};
    boxL[] += {newl}; Circle(boxL[11]) = {boxP[4],boxP[0],boxP[6]};

    ll = newll; Line Loop(ll) = {boxL[0],boxL[6],-boxL[4]};
    boxS[] += {news}; Ruled Surface(boxS[0]) = {ll};
    ll = newll; Line Loop(ll) = {boxL[1],boxL[8],-boxL[6]};
    boxS[] += {news}; Ruled Surface(boxS[1]) = {ll};
    ll = newll; Line Loop(ll) = {boxL[2],boxL[10],-boxL[8]};
    boxS[] += {news}; Ruled Surface(boxS[2]) = {ll};
    ll = newll; Line Loop(ll) = {boxL[3],boxL[4],-boxL[10]};
    boxS[] += {news}; Ruled Surface(boxS[3]) = {ll};
    ll = newll; Line Loop(ll) = {boxL[0],boxL[7],-boxL[5]};
    boxS[] += {news}; Ruled Surface(boxS[4]) = {ll};
    ll = newll; Line Loop(ll) = {boxL[1],boxL[9],-boxL[7]};
    boxS[] += {news}; Ruled Surface(boxS[5]) = {ll};
    ll = newll; Line Loop(ll) = {boxL[2],boxL[11],-boxL[9]};
    boxS[] += {news}; Ruled Surface(boxS[6]) = {ll};
    ll = newll; Line Loop(ll) = {boxL[3],boxL[5],-boxL[11]};
    boxS[] += {news}; Ruled Surface(boxS[7]) = {ll};

    sl = newsl; Surface Loop(sl) = {boxS[]};
    boxV = newv; Volume(boxV) = {sl};
Return

//////////////////// END MACRO //////////////////////

Ri=RI;
Re=RE;
Z=ZB;
L=Zs;
Call cyl;
baseP[] = cylP[];
baseL[] = cylL[];
baseS[] = cylS[];
baseV[] = cylV[];

Z=Ze;
L=ZT-Z;
Call cyl;
topP[] = cylP[];
topL[] = cylL[];
topS[] = cylS[];
topV[] = cylV[];

baseIP = baseP[2];
baseEP = baseP[3];
topIP = topP[0];
topEP = topP[1];
baseLine = baseL[1];
topLine = topL[0];
Call Helix;
hel1P[] = helP[];
hel1L[] = helL[];
hel1S[] = helS[];
hel1V[] = helV[];

Ri=-RI;
Re=-RE;
baseIP = baseP[4*(n/2)+2];
baseEP = baseP[4*(n/2)+3];
topIP = topP[4*(n/2)+0];
topEP = topP[4*(n/2)+1];
baseLine = baseL[8*(n/2)+1];
topLine = topL[8*(n/2)];
Call Helix;
hel2P[] = helP[];
hel2L[] = helL[];
hel2S[] = helS[];
hel2V[] = helV[];

Call CondDown;
Call CondUp;
Call CondBet;

Call Box;

For i In { 1 : n }
    base[] += {baseS[5*i]};
    top[] += {topS[5*i+1]};
    interior[] += {baseS[5*i+2],topS[5*i+2]};
    exterior[] += {baseS[5*i+3],topS[5*i+3]};
    cond[] += {baseV[i-1],topV[i-1]};
EndFor
For i In { 0 : n/2 }
    interior[] += {condDS[6*i+2],condDS[6*i+3],condUS[6*i+2],condUS[6*i+3]};
    exterior[] += {condDS[6*i+4],condDS[6*i+5],condUS[6*i+4],condUS[6*i+5]};
    cond[] += {condDV[2*i],condDV[2*i+1],condUV[2*i],condUV[2*i+1]};
EndFor
For i In { 1 : n*nbTour-n/2 }
    interior[] += {condBS[6*i+2],condBS[6*i+3]};
    exterior[] += {condBS[6*i+4],condBS[6*i+5]};
    cond[] += {condBV[2*(i-1)],condBV[2*(i-1)+1]};
EndFor
For i In { 1 : n*nbTour }
    helixS1[] += {hel1S[5*i+1],hel1S[5*i+2]};
    helixS2[] += {hel2S[5*i+1],hel2S[5*i+2]};
    helixV1[] += {hel1V[i-1]};
    helixV2[] += {hel2V[i-1]};
EndFor
helixS1[] += {hel1S[0],hel1S[5*n*nbTour]};
helixS2[] += {hel2S[0],hel2S[5*n*nbTour]};

Physical Volume("cond") = {cond[]};
Physical Volume("helix1") = {helixV1[]};
Physical Volume("helix2") = {helixV2[]};
Physical Volume("box") =  {boxV};
Physical Surface("base") = {base[]};
Physical Surface("top") = {top[]};
Physical Surface("int") = {interior[]};
Physical Surface("ext") = {exterior[]};
Physical Surface("helixS1") = {helixS1[]};
Physical Surface("helixS2") = {helixS2[]};
