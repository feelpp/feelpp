#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <map>

#include <qt3/qfile.h>
#include <qt3/qregexp.h>
#include <qt3/qdir.h>

#include <adfstream.hpp>

void
generate_binary_functions( std::string const& __prefix )
{
    struct
    {
        const char* funcObj;
        const char* funcName;
        const char* comment;
        const char* grad;
        const char* hess;
        const char* grad_std;
        const char* hess_std;
    } funcs [] =
    {
        {
            "AdFuncPow", "pow", "Power function",
            "expr2_.value() * expr1_.grad(__i) * pow(expr1_.value(),expr2_.value()-1)",
            "expr2_.value() * ( (expr2_.value()-1) * expr1_.grad(__i) * expr1_.grad(__j)  * pow(expr1_.value(),expr2_.value()-2) +  expr1_.hessian(__i,__j) * pow(expr1_.value(),expr2_.value()-1))",
            "expr2_.value() * expr1_.grad(__i) * std::pow(expr1_.value(),expr2_.value()-1)",
            "expr2_.value() * ( (expr2_.value()-1) * expr1_.grad(__i) * expr1_.grad(__j) * std::pow(expr1_.value(),expr2_.value()-2) + expr1_.hessian(__i,__j) * std::pow(expr1_.value(),expr2_.value()-1))"
        },
    };

    std::cout << "Generating <adbinaryfunctions.hpp>" << std::endl;

    //OperandTuple operands(2);

    Feel::ADOfstream ofs( "adbinaryfunctions.hpp", "ADType expression templates", __FILE__, "AD_BINARY_FUNCS_HPP" );
    ofs.beginNamespace( "Feel" );
    const int numOperators = 1;   // Should be 18

    for ( int i=0; i < numOperators; ++i )
    {
        ofs << "/****************************************************************************" << std::endl
            << " * " << funcs[i].comment << std::endl
            << " ****************************************************************************/" << std::endl;

        QStringList lines;
        std::string fname = __prefix + "/adbinaryfunctions.tmpl.hpp";
        QFile file( fname.c_str() );

        if ( file.exists() && file.open( IO_ReadOnly ) )
        {
            QTextStream stream( &file );
            QString line;

            while ( !stream.eof() )
            {
                line = stream.readLine(); // line of text excluding '\n'

                QRegExp rx;
                rx.setPattern ( "@NAME@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].funcObj;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@FCT@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].funcName;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@GRDI@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].grad;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@HESSIJ@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].hess;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@GRDI_STD@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].grad_std;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@HESSIJ_STD@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].hess_std;
                    line.replace( rx, __os.str().c_str() );
                }

                ofs << line << std::endl;
            }

            file.close();
        }
    }

    ofs.endNamespace();
}

void
generate_functions( std::string const& __prefix )
{
    struct
    {
        const char* funcObj;
        const char* funcName;
        const char* comment;
        const char* grad;
        const char* hess;
        const char* grad_std;
        const char* hess_std;
    } funcs [] =
    {
        {
            "ADFuncCos", "cos",  "Cosinus function",
            "-expr_.grad(__i)*sin( expr_.value() )",
            "-expr_.hessian(__i,__j)*sin( expr_.value() ) - expr_.grad( __i )*expr_.grad( __j ) * cos( expr_.value() )",
            "-expr_.grad(__i)*std::sin( expr_.value() )",
            "-expr_.hessian(__i,__j)*std::sin( expr_.value() ) - expr_.grad( __i )*expr_.grad( __j ) * std::cos( expr_.value() )"
        },
        {
            "ADFuncSin", "sin", "Sinus function",
            "expr_.grad( __i )*cos(expr_.value())",
            "expr_.hessian(__i,__j)*cos( expr_.value() ) - expr_.grad( __i )*expr_.grad( __j ) * sin( expr_.value() )",
            "expr_.grad( __i )*std::cos(expr_.value())",
            "expr_.hessian(__i,__j)*std::cos( expr_.value() ) - expr_.grad( __i )*expr_.grad( __j ) * std::sin( expr_.value() )",
        },
        {
            "ADFuncTan", "tan", "Tangent function",
            "expr_.grad(__i)*(1.+tan(expr_.value())*tan(expr_.value()))",
            "2*tan(expr_.value())*(1+tan(expr_.value())*tan(expr_.value()))*expr_.grad(__i)*expr_.grad(__j) + (1+tan(expr_.value())*tan(expr_.value()))*expr_.hessian(__i,__j)",
            "expr_.grad(__i)*(1.+std::tan(expr_.value())*std::tan(expr_.value()))",
            "2*std::tan(expr_.value())*(1+std::tan(expr_.value())*std::tan(expr_.value()))*expr_.grad(__i)*expr_.grad(__j) + (1+std::tan(expr_.value())*std::tan(expr_.value()))*expr_.hessian(__i,__j)"
        },
        {
            "ADFuncSqrt", "sqrt", "Sqrt function",
            "expr_.grad(__i)/(2.*sqrt(expr_.value()))",
            "expr_.hessian(__i,__j)/(2.*sqrt(expr_.value()))-expr_.grad(__i)*expr_.grad(__j)/(4.*pow(expr_.value(),1.5) )",
            "expr_.grad(__i)/(2.*std::sqrt(expr_.value()))",
            "expr_.hessian(__i,__j)/(2.*std::sqrt(expr_.value()))-expr_.grad(__i)*expr_.grad(__j)/(4.*std::pow(expr_.value(),1.5) )"
        },
        {
            "AdFuncExp", "exp", "Exponential function",
            "expr_.grad(__i)*exp(expr_.value())",
            "expr_.hessian(__i,__j)*exp(expr_.value()) + expr_.grad(__i)*expr_.grad(__j)*exp(expr_.value())",
            "expr_.grad(__i)*std::exp(expr_.value())",
            "expr_.hessian(__i,__j)*std::exp(expr_.value()) + expr_.grad(__i)*expr_.grad(__j)*std::exp(expr_.value())"
        },
        {
            "AdFuncLog", "log", "Logarithm function",
            "expr_.grad(__i)/expr_.value()",
            "expr_.hessian(__i,__j)/expr_.value() - expr_.grad(__i)*expr_.grad(__j)/(pow(expr_.value(),2.))",
            "expr_.grad(__i)/expr_.value()",
            "expr_.hessian(__i,__j)/expr_.value() - expr_.grad(__i)*expr_.grad(__j)/(std::pow(expr_.value(),2.))"
        },
        {
            "AdFuncLog10", "log10", "Log10 function",
            "expr_.grad(__i)/(log(value_type(10))*expr_.value())",
            "expr_.hessian(__i,__j)/(log(value_type(10))*expr_.value()) -  expr_.grad(__i)*expr_.grad(__j)/(log(value_type(10))*pow(expr_.value(),2.))",
            "expr_.grad(__i)/(std::log(value_type(10))*expr_.value())",
            "expr_.hessian(__i,__j)/(std::log(value_type(10))*expr_.value()) -  expr_.grad(__i)*expr_.grad(__j)/(std::log(value_type(10))*std::pow(expr_.value(),2.))"
        },
        {
            "AdFuncAbs", "abs", "Absolute function",
            "signum(expr_.value()) * expr_.grad(__i)",
            "signum(expr_.value()) * expr_.hessian(__i,__j)",
            "signum(expr_.value()) * expr_.grad(__i)",
            "signum(expr_.value()) * expr_.hessian(__i,__j)"
        },
    };

    std::cout << "Generating <adfunctions.hpp>" << std::endl;

    //OperandTuple operands(2);

    Feel::ADOfstream ofs( "adfunctions.hpp", "ADType expression templates",
                          __FILE__, "AD_FUNCS_HPP" );
    ofs.beginNamespace( "Feel" );
    const int numOperators = 8;   // Should be 18

    for ( int i=0; i < numOperators; ++i )
    {
        ofs << "/****************************************************************************" << std::endl
            << " * " << funcs[i].comment << std::endl
            << " ****************************************************************************/" << std::endl;

        QStringList lines;
        std::string fname = __prefix + "/adfunctions.tmpl.hpp";
        QFile file( fname.c_str() );

        if ( file.exists() && file.open( IO_ReadOnly ) )
        {
            QTextStream stream( &file );
            QString line;

            while ( !stream.eof() )
            {
                line = stream.readLine(); // line of text excluding '\n'

                QRegExp rx;
                rx.setPattern ( "@NAME@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].funcObj;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@FCT@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].funcName;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@GRDI@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].grad;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@HESSIJ@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].hess;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@GRDI_STD@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].grad_std;
                    line.replace( rx, __os.str().c_str() );
                }

                rx.setPattern ( "@HESSIJ_STD@" );

                if ( rx.search( line ) != -1 )
                {
                    std::ostringstream __os;
                    __os << funcs[i].hess_std;
                    line.replace( rx, __os.str().c_str() );
                }

                ofs << line << std::endl;
            }

            file.close();
        }
    }

    ofs.endNamespace();
}
void generate_binary_operators( std::string const& __prefix )
{

    std::cout << "Generating <adbinaryoperators.hpp>" << std::endl;

    //OperandTuple operands(2);

    Feel::ADOfstream ofs( "adbinaryoperators.hpp", "ADType expression templates (2 operands)",
                          __FILE__, "AD_BOPS_HPP" );

    struct
    {
        const char* opSymbol;
        bool        nonIntOperands;
        bool        complexOperands;
        const char* opApplicName;
        const char* comment;
    } ops[] =
    {
        { "+",  true,  true,  "Add",            "Addition Operators" },
        { "-",  true,  true,  "Subtract",       "Subtraction Operators" },
        { "*",  true,  true,  "Multiply",       "Multiplication Operators" },
        { "/",  true,  true,  "Divide",         "Division Operators" },
        { "^",  true,  true,  "Pow",            "Pow Operators" },
        { ">",  true,  true,  "Greater",        "Greater-than Operators" },
        { "<",  true,  false, "Less",           "Less-than Operators" },
        { ">=", true,  false, "GreaterOrEqual", "Greater or equal (>=) operators" },
        { "<=", true,  false, "LessOrEqual",    "Less or equal (<=) operators" },
        { "==", true,  true,  "Equal",          "Equality operators" },
        { "!=", true,  true,  "NotEqual",       "Not-equal operators" },
    };

#if 0
    { "%",  false, false, "Modulo",         "Modulus Operators" },


    { "==", true,  true,  "Equal",          "Equality operators" },
    { "!=", true,  true,  "NotEqual",       "Not-equal operators" },
    { "&&", false, false, "LogicalAnd",     "Logical AND operators" },
    { "||", false, false, "LogicalOr",      "Logical OR operators" },



    { "&",  false, false, "BitwiseAnd",     "Bitwise And Operators" },
    { "|",  false, false, "BitwiseOr",      "Bitwise Or Operators" },
    { ">>", false, false, "ShiftRight",     "Shift right Operators" },
    { "<<", false, false, "ShiftLeft",      "Shift left Operators" }
};
#endif
const int numOperators = 11;   // Should be 18

for ( int i=0; i < numOperators; ++i )
{
    ofs << "/****************************************************************************" << std::endl
        << " * " << ops[i].comment << std::endl
        << " ****************************************************************************/" << std::endl;
    std::string __fname = __prefix + "/adbinaryoperators.tmpl.hpp";
    std::cout << "loading template " << __fname << "\n";
    QStringList lines;
    QFile file( __fname.c_str() );

    if ( file.exists() && file.open( IO_ReadOnly ) )
    {
        QTextStream stream( &file );
        QString line;

        while ( !stream.eof() )
        {
            line = stream.readLine(); // line of text excluding '\n'

            QRegExp rx;
            rx.setPattern ( "@OP@" );

            if ( rx.search( line ) != -1 )
            {
                std::ostringstream __os;
                __os << "operator" << ops[i].opSymbol;
                line.replace( rx, __os.str().c_str() );
            }

            rx.setPattern ( "@STYPE@" );

            if ( rx.search( line ) != -1 )
            {
                std::ostringstream __os;
                __os << "ADBinary" << ops[i].opApplicName;
                line.replace( rx, __os.str().c_str() );
            }

            ofs << line << std::endl;
        }

        file.close();
    }
}


}

int main( int argc, char** argv )
{
    generate_binary_operators( argv[1] );
    generate_functions( argv[1] );
    generate_binary_functions( argv[1] );
}

