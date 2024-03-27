#include <thread>
#include <vector>
#include <array>
#include <typeinfo>
#include <iostream>
#include <mutex>
#include <sched.h>
#include <pthread.h>

#include<algorithm> //for Each_fors



#include <string>

#include <utility>
#include <functional>

#include <future>

#include <cassert>

#include <chrono>

#include <type_traits>

#include <list>
#include <ranges>


//================================================================================================================================

#define clrscr() printf("\033[H\033[2J")
#define color(param) printf("\033[%sm",param)


//================================================================================================================================
// Tools to manage the console: colors, cursor,...

/* 
    0  réinitialisation         1  haute intensité (des caractères)
    5  clignotement             7  video inversé
    30, 31, 32, 33, 34, 35, 36, 37 couleur des caractères
    40, 41, 42, 43, 44, 45, 46, 47 couleur du fond
    noir, rouge, vert, jaune, bleu, magenta, cyan et blanc 
    0,      1,    2,     3,     4,     5,      6,      7
    color("40;37")               
*/


void Color(int number)
{
    if (number>7) { number=number % 8; }
    number+=30;
    printf("\033[%im",number);
}

void ColorBackground(int number)
{
    if (number>7) { number=number % 8; }
    number+=40;
    printf("\033[%im",number);
}

//BEGIN::CONSOLE TOOLS FOR LINUX
void CONSOLE_ClearScreen()                     { printf("\033[2J"); }
void CONSOLE_SaveCursorPosition()              { printf("\033[s"); }
void CONSOLE_RestoreCursorPosition()           { printf("\033[u");}
void CONSOLE_SetCursorPosition(int x, int y)   { printf("\033[%d;%dH",y+1,x+1); }
void CONSOLE_GetCursorPosition(int* x, int* y) { printf("\033[6n");  int err=scanf("\033[%d;%dR", x, y); }
void CONSOLE_CursorBlinkingOnOff()             { printf("\033[?12l"); }
void CONSOLE_CursorHidden(bool q)              { q ? printf("\e[?25l"):printf("\e[?25h"); }

void CONSOLE_PRINT_PERCENTAGE(int x, int y, int k, int nb)
{
	float v = ((float(k)) / float(nb)) * 100.0;
	CONSOLE_SetCursorPosition(x,y); printf(" %3.1f %%", v);
}

void CONSOLE_PRINT_PERCENTAGE(int k, int nb)
{
	float v = ((float(k)) / float(nb)) * 100.0;
	CONSOLE_RestoreCursorPosition(); printf(" %3.1f %%", v);
}

void CONSOLE_Color(int n)
{
    switch(n) {
        case 0:
            printf("\e]P0000000"); ///black
        break;
        case 1:
            printf("\e]P1D75F5F"); //darkred
        break;
        case 2:
            printf("\e]P287AF5F"); //darkgreen
        break;
        case 3:
            printf("\e]P3D7AF87"); //brown
        break;
        case 4:
            printf("\e]P48787AF"); //darkblue
        break;
        case 5:
            printf("\e]P5BD53A5"); //darkmagenta
        break;
        case 6:
            printf("\e]P65FAFAF"); //darkcyan
        break;
        case 7:
            printf("\e]P7E5E5E5"); //lightgrey
        break;
        case 8:
            printf("\e]P82B2B2B"); //darkgrey
        break;
        case 9:
            printf("\e]P9E33636"); //red
        break;
        case 10:
            printf("\e]PA98E34D"); //green
        break;
        case 11:
            printf("\e]PBFFD75F"); //yellow
        break;
        case 12:
            printf("\e]PC7373C9"); //blue
        break;
        case 13:
            printf("\e]PDD633B2"); //magenta
        break;
        case 14:
            printf("\e]PE44C9C9"); //cyan
        break;
        case 15:
            printf("\e]PFFFFFFF"); //white
        break;
    }
}
//1-CONSOLE_SaveCursorPosition();
//2-CONSOLE_CursorHidden(true); 
//3-LOOP   CONSOLE_PRINT_PERCENTAGE(i+1,M_listFaceMarkers.size());
//4-CONSOLE_CursorHidden(false); 
//END::CONSOLE TOOLS



//================================================================================================================================
// Functions that provide the characteristics, type of a variable,...

template <typename T>
std::string GetCallTyp(T& fcv)
{
    std::string rep="Unknow";
    std::string s1=typeid(fcv).name();
    bool qconst=std::is_const_v<T>;
    bool qclass=std::is_class_v<T>;
    bool qfunction=std::is_function_v<T>;
    bool qpointer=std::is_pointer_v<T>;
    bool qarray=std::is_array_v<T>;
    bool qenum=std::is_enum_v<T>;
    bool qunion=std::is_union_v<T>;
    bool qvoid=std::is_void_v<T>;

    std::cout << "type : " << s1 <<'\n';
    std::cout << "void : " << qvoid <<" ";  
    std::cout << "enum : " << qenum <<" ";
    std::cout << "union : " << qunion <<" "; 
    std::cout << "const : " << qconst <<" ";

    std::cout << "array : " << qarray <<" ";
    std::cout << "pointer : " << qpointer <<" ";
    std::cout << "function : " << qfunction <<" ";
    std::cout << "class : " << qclass << '\n';

    int l=s1.length();
    if (l>1) {
        if (s1.find("16SpScalarDataModeIL16SpDataAccessMode0") != std::string::npos) { rep="SpRead"; }
        if (s1.find("16SpScalarDataModeIL16SpDataAccessMode1") != std::string::npos) { rep="SpWrite"; }
        if (s1.find("16SpScalarDataModeIL16SpDataAccessMode3") != std::string::npos) { rep="SpCommutativeWrite"; }
        if (s1.find("17SpCallableWrapperILb1ERKiL14SpCallableType0EE") != std::string::npos) { rep="SpCpu"; }
        if (s1.find("17SpCallableWrapperILb0ERKiL14SpCallableType1EE") != std::string::npos) { rep="SpCuda"; }
        if (s1.find("17SpCallableWrapperILb0ERKiL14SpCallableType2EE") != std::string::npos) { rep="SpHip"; }
        if (s1.find("vector") != std::string::npos) { 
            if (qconst) { rep="Const Vector"; } else { rep="Vector"; }
        }
    }
    else 
    {
        if (s1.find("i")!= std::string::npos) { 
            if (qconst) { rep="Const Integer"; } else { rep="Integer"; }
        }
        if (s1.find("f")!= std::string::npos) { 
            if (qconst) { rep="Const Float"; } else { rep="Float"; }
        }
        if (s1.find("d")!= std::string::npos) { 
            if (qconst) { rep="Const Double"; } else { rep="Double "; }
         }
    }
    return (rep);
}


template <typename T>
int GetSignatureSpecxFunction(T& fcv)
{
    int res=0;
    std::string rep="Unknow";
    std::string s1=typeid(fcv).name();
    bool qconst=std::is_const_v<T>;
    bool qclass=std::is_class_v<T>;
    bool qfunction=std::is_function_v<T>;
    bool qpointer=std::is_pointer_v<T>;
    bool qarray=std::is_array_v<T>;
    bool qenum=std::is_enum_v<T>;
    bool qunion=std::is_union_v<T>;
    bool qvoid=std::is_void_v<T>;

    int l=s1.length();
    if (l>1) {
        if (s1.find("16SpScalarDataModeIL16SpDataAccessMode0") != std::string::npos) { res=1; } //SpRead
        if (s1.find("16SpScalarDataModeIL16SpDataAccessMode1") != std::string::npos) { res=2; } //SpWrite
        if (s1.find("16SpScalarDataModeIL16SpDataAccessMode3") != std::string::npos) { res=3; } //SpCommutativeWrite
    }
    return (res);
}


//================================================================================================================================

auto GetListNameObjects(std::string ChName)
{
    std::ifstream FICH(ChName);
    std::string lineA;
    std::vector<std::string> ListObjects;
    while (getline(FICH, lineA)) { ListObjects.push_back(lineA); } 
    FICH.close();
    return ListObjects;
}

//================================================================================================================================

void WriteMat(int Mat[][3],int row,int col)
{
    for(int i = 0; i < row; i++) {
        for(int j = 0; j < col; j++) {
            Color(j+1);
            std::cout << Mat[i][j];
            std::cout << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    Color(7);
}

//================================================================================================================================
