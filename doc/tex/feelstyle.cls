\def\filedate{2003/10/29}
\def\fileversion{1}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{lifestyle}[2003/10/29]
\message{`lifestyle v\fileversion, \filedate\space (Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>)}

\RequirePackage{ifthen}

\newboolean{life@article} \setboolean{life@article}{true}
\newboolean{life@report}  \setboolean{life@report}{false}
\newboolean{life@toc}     \setboolean{life@toc}{false}
\newboolean{life@stretch} \setboolean{life@stretch}{false}

\DeclareOption{article}  {\setboolean{life@article}{true}\setboolean{life@report}{false}}
\DeclareOption{report}   {\setboolean{life@report}{true}\setboolean{life@article}{false}}
\DeclareOption{toc}      {\setboolean{life@toc}{true}}
\DeclareOption{stretch}  {\setboolean{life@stretch}{true}}

\ifthenelse{\boolean{life@article}}
{
  \DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
}{
  \ifthenelse{\boolean{life@report}}
  {
    \DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
  }{
    %% default is article
    \DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
  }
}

%%
%% Execution of Options
%%

\ProcessOptions

\ifthenelse{\boolean{life@article}}
{
  \LoadClass[titlepage]{article}
}{
  \ifthenelse{\boolean{life@report}}
  {
    \LoadClass[titlepage]{report}
  }{
    %% default is article
    \LoadClass[titlepage]{article}
  }
}



\RequirePackage{a4wide}
\RequirePackage{url}
\RequirePackage{fancyhdr}
\RequirePackage{fancybox}
\RequirePackage{xspace}
\RequirePackage{multicol}
\RequirePackage{makeidx}
\RequirePackage{subfigure}
\RequirePackage{verbatim}
\RequirePackage{rcs}
\RequirePackage{colortbl}

\RequirePackage[english]{babel}
\RequirePackage{tabularx}

%\RequirePackage[pdftex,colorlinks]{hyperref}

% Graphics rules
%
\RequirePackage{color}
\RequirePackage{graphicx}
\DeclareGraphicsExtensions{.jpg,.png,.pdf,.ps}
\DeclareGraphicsRule{*}{mps}{*}{}

% Some maths
%
\RequirePackage{amsmath,amsfonts,amssymb,amscd,amsthm,amsbsy}

%%\newcommand{\lifedraftdate}{\small Revision: \RCSRevision, \RCSDate, \RCSTime\ UTC Printed: \today}

\newenvironment{rationale}{\renewcommand{\labelitemi}{\textit{Rationale:}}
  \begin{itemize}
  \item
    }{  \end{itemize}}

%%
%% Index
%%
\newcommand{\ix}[1]{#1\index{#1}}
\newcommand{\ixn}[1]{\index{#1}}
\newcommand{\ixv}[1]{{\tt #1}\index{#1}}
\newcommand{\ixs}[2]{#1\index{#1|see{#2}}\index{#2!#1}}
\newcommand{\ixns}[2]{\index{#1|see{#2}}\index{#2!#1}}
\newcommand{\ixvs}[2]{{\tt #1}\index{#1|see{#2}}\index{#2!#1}}
\newcommand{\nindex}[1]{\textit{#1}\index[name]{#1}}


\newcommand{\cppobj}[1]{\index{C++ Objects!#1}}
\newcommand{\lifetitle}[2]{\raggedright{\Huge #1}\\
 \rule[2ex]{\linewidth}{1mm}\\
 {\hfill\lifedraftdate}\\
 \vfill{\Large #2}\\
 \rule[2ex]{\linewidth}{.5mm}}

\newcommand{\thelibrary}{\textsl{LIFE}\xspace}
\newcommand{\life}{\textsl{Life}\xspace}

%%
%% generate hints
%%
\usepackage{fancybox,framed}
\newcommand{\hinttitle}[1]{\mbox{#1}\xspace}
\definecolor{hintcolor}{rgb}{.95,1,.95}
% \newenvironment{\hint}[1]{
% \begin{framed}
% \colorbox{hintcolor}{\parbox{.65\linewidth}{
% % \begin{tabular}[.65\linewidth]{l}
% % \hinttitle{#1}\\
% % \rule[2ex]{\linewidth}{1mm}\\
% }{\end{framed}}
\makeatletter
\newenvironment{hint}[1]{%

  \def\FrameCommand{%
    % \setlength{\fboxrule}{\FrameRule}%
    % \setlength{\fboxsep}{\FrameSep}%
    \shadowbox
  }%
  \fram@d{\advance\hsize-\width}%
  #1 Hint\\
  \rule[2ex]{\linewidth}{.5mm}\\
}{%
  \endfram@d
}
\makeatother




\makeglossary
\makeindex

\AtBeginDocument{

  \maketitle

  \ifthenelse{\boolean{life@stretch}}{\baselineskip 0.33in}{}

  \ifthenelse{\boolean{life@toc}}{
    \setlength{\columnseprule}{0.6pt}
    \begin{multicols}{2}
      \tableofcontents
    \end{multicols}
  }

}