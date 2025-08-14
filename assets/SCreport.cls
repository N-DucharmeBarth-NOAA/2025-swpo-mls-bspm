\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{SCreport}[2025/03/17 v20.4 WCPFC Scientific Committee report]

% ----------
% 1  Inherit
% ----------

% 1.1 Class
\LoadClass{article}

% ---------------
% 2  Dependencies
% ---------------

% 2.1 Packages
\RequirePackage{setspace}
\RequirePackage{soul}
\RequirePackage{tocloft}
\RequirePackage{fvextra}
\RequirePackage[usenames,dvipsnames,svgnames]{xcolor}
\RequirePackage{xspace}
\RequirePackage{hyperref}


% -----------
% 3  Settings
% -----------

% 3.1 Color
\definecolor{SCblue}{rgb}{0,0.1,0.7}
\hypersetup{colorlinks,linkcolor=SCblue,filecolor=SCblue,urlcolor=SCblue,citecolor=SCblue}
\def\hldefaultcol{Khaki}
\sethlcolor{\hldefaultcol}

% 3.2 Spacing
\frenchspacing
\onehalfspacing
\setlength\cftbeforesecskip{0pt}  % vertical space in TOC

% 3.3 Encoding
\DeclareUnicodeCharacter{2212}{-}  % convert minus sign to hyphen

% 3.4 Inline code formatting
\newcommand{\inlinecode}[1]{\colorbox[HTML]{f6f8fa}{\ttfamily\textcolor[HTML]{6f42c1}{\detokenize{#1}}}}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{commandchars=\\\{\},fontsize=\small,frame=single,rulecolor=\color{gray!30},framesep=2mm,baselinestretch=1.2,bgcolor=\color[HTML]{f6f8fa}}

% -------------------------
% 4  Commands and variables
% -------------------------

% 4.1 Misc
\newcommand{\LnS}{\mathrm{log}\hspace{0.4mm}}
\newcommand{\degree}{$^{\circ}$\xspace}
\newcommand{\ra}{$\rightarrow{\;}$}
\newcommand{\vs}{\emph{vs.}\xspace}
\newcommand{\via}{\emph{via}\xspace}
\newcommand{\done}[1]{\st{#1}} % strike text

% 4.2 SPC
\newcommand{\spc}{Oceanic Fisheries Programme of the Pacific Community}
\newcommand{\lrp}{$20\%\mbox{\emph{SB}}_{F=0}$\xspace}
\newcommand{\mfcl}{MULTIFAN-CL\xspace}
\newcommand{\sss}{Stock Synthesis\xspace}
\newcommand{\spalb}{south Pacific albacore tuna\xspace}
\newcommand{\alb}{albacore tuna\xspace}
\newcommand{\bet}{bigeye tuna\xspace}
\newcommand{\swpo}{Southwest Pacific Ocean\xspace}
\newcommand{\mls}{striped marlin\xspace}

% 4.3 Spawning biomass
\newcommand{\SB}{$\mbox{\emph{SB}}$\xspace}
\newcommand{\sbcurr}{$\mbox{\emph{SB}}_\mathrm{current}$\xspace}
\newcommand{\sbcsbfo}{$\mbox{\emph{SB}}_\mathrm{current}/\mbox{\emph{SB}}_{F=0}$\xspace}
\newcommand{\sbcsbmsy}{$\mbox{\emph{SB}}_\mathrm{current}/\mbox{\emph{SB}}_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\sbfo}{$\mbox{\emph{SB}}_{F=0}$\xspace}
\newcommand{\sblatest}{$\mbox{\emph{SB}}_\mathrm{latest}$\xspace}
\newcommand{\sblsbfo}{$\mbox{\emph{SB}}_\mathrm{latest}/\mbox{\emph{SB}}_{F=0}$\xspace}
\newcommand{\sblsbmsy}{$\mbox{\emph{SB}}_\mathrm{latest}/\mbox{\emph{SB}}_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\sblsbo}{$\mbox{\emph{SB}}_\mathrm{latest}/\mbox{\emph{SB}}_{0}$\xspace}
\newcommand{\sbmsy}{$\mbox{\emph{SB}}_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\sbmsysbo}{$\mbox{\emph{SB}}_{\mbox{\smaller[3]{MSY}}}/\mbox{\emph{SB}}_{0}$\xspace}
\newcommand{\sbmsysbfo}{$\mbox{\emph{SB}}_{\mbox{\smaller[3]{MSY}}}/\mbox{\emph{SB}}_{F=0}$\xspace}
\newcommand{\sbo}{$\mbox{\emph{SB}}_{0}$\xspace}
\newcommand{\sbrecent}{$\mbox{\emph{SB}}_\mathrm{recent}$\xspace}
\newcommand{\sbrsbfo}{$\mbox{\emph{SB}}_\mathrm{recent}/\mbox{\emph{SB}}_{F=0}$\xspace}
\newcommand{\sbrsbmsy}{$\mbox{\emph{SB}}_\mathrm{recent}/\mbox{\emph{SB}}_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\sbsbfo}{$\mbox{\emph{SB}}/\mbox{\emph{SB}}_{F=0}$\xspace}
\newcommand{\sbsbmsy}{$\mbox{\emph{SB}}/\mbox{\emph{SB}}_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\sbt}{$\mbox{\emph{SB}}_t$\xspace}
\newcommand{\sbtfo}{$\mbox{\emph{SB}}_{t,F=0}$\xspace}
\newcommand{\sbtsbfo}{$\mbox{\emph{SB}}_t/\mbox{\emph{SB}}_{t,F=0}$\xspace}
\newcommand{\sbtsbmsy}{$\mbox{\emph{SB}}_t/\mbox{\emph{SB}}_{\mbox{\smaller[3]{MSY}}}$\xspace}

% 4.4 Biomass
\newcommand{\bcurr}{$\mbox{\emph{B}}_\mathrm{current}$\xspace}
\newcommand{\bcbfo}{$\mbox{\emph{B}}_\mathrm{current}/\mbox{\emph{B}}_{F=0}$\xspace}
\newcommand{\bcbmsy}{$\mbox{\emph{B}}_\mathrm{current}/\mbox{\emph{B}}_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\bfo}{$\mbox{\emph{B}}_{F=0}$\xspace}
\newcommand{\blatest}{$\mbox{\emph{B}}_\mathrm{latest}$\xspace}
\newcommand{\blbfo}{$\mbox{\emph{B}}_\mathrm{latest}/\mbox{\emph{B}}_{F=0}$\xspace}
\newcommand{\blbmsy}{$\mbox{\emph{B}}_\mathrm{latest}/\mbox{\emph{B}}_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\bmsy}{$\mbox{\emph{B}}_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\bo}{$\mbox{\emph{B}}_{0}$\xspace}

% 4.5 Fishing mortality and catch
\newcommand{\clatest}{$\mbox{\emph{C}}_\mathrm{latest}$\xspace}
\newcommand{\clatestmsy}{$\mbox{\emph{C}}_\mathrm{latest}/F_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\fcurr}{$F_\mathrm{current}$\xspace}
\newcommand{\fmsy}{$F_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\ffmsy}{$F/F_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\fmult}{$f_\mathrm{mult}$\xspace}
\newcommand{\frecent}{$F_\mathrm{recent}$\xspace}
\newcommand{\fref}{$F_\mathrm{recent}/F_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\ftfmsy}{$F_t/F_{\mbox{\smaller[3]{MSY}}}$\xspace}
\newcommand{\msy}{$\mbox{MSY}$\xspace}
\newcommand{\yfcurr}{$Y_{F_\mathrm{current}}$\xspace}
\newcommand{\yfrecent}{$Y_{F_\mathrm{recent}}$\xspace}

% 4.6 Colors
\newcommand{\hllb}[1]{\sethlcolor{PaleTurquoise} \hl{#1} \sethlcolor{\hldefaultcol}}
\newcommand{\hlam}[1]{\sethlcolor{Aquamarine} \hl{#1} \sethlcolor{\hldefaultcol}}
\newcommand{\hlcyan}[1]{\sethlcolor{Cyan} \hl{#1} \sethlcolor{\hldefaultcol}}
\newcommand{\hlyel}[1]{\sethlcolor{Gold} \hl{#1} \sethlcolor{\hldefaultcol}}
\newcommand{\hlred}[1]{\sethlcolor{OrangeRed} \hl{#1} \sethlcolor{\hldefaultcol}}
\newcommand{\hlgreen}[1]{\sethlcolor{YellowGreen} \hl{#1} \sethlcolor{\hldefaultcol}}
\newcommand{\hlcoral}[1]{\sethlcolor{Coral} \hl{#1} \sethlcolor{\hldefaultcol}}
