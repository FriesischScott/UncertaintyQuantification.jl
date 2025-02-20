import{_ as e,c as t,o as i,j as s,ai as n,a as l}from"./chunks/framework.BuyV_eSC.js";const h="/UncertaintyQuantification.jl/v0.12.0/assets/stiffness-samples.BanvBlFC.svg",p="/UncertaintyQuantification.jl/v0.12.0/assets/stiffness-point-estimate-uniform.B0P0M9b-.svg",T="/UncertaintyQuantification.jl/v0.12.0/assets/stiffness-point-estimate-normal.DAWG9op6.svg",a1=JSON.parse('{"title":"Bayesian Updating","description":"","frontmatter":{},"headers":[],"relativePath":"examples/bayesianupdating.md","filePath":"examples/bayesianupdating.md","lastUpdated":null}'),d={name:"examples/bayesianupdating.md"},r={class:"MathJax",jax:"SVG",display:"true",style:{direction:"ltr",display:"block","text-align":"center",margin:"1em 0",position:"relative"}},Q={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-2.149ex"},xmlns:"http://www.w3.org/2000/svg",width:"15.324ex",height:"5.43ex",role:"img",focusable:"false",viewBox:"0 -1450 6773.1 2400","aria-hidden":"true"},k={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},o={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"7.957ex",height:"1.934ex",role:"img",focusable:"false",viewBox:"0 -705 3517.1 855","aria-hidden":"true"},m={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},g={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"7.957ex",height:"1.934ex",role:"img",focusable:"false",viewBox:"0 -705 3517.1 855","aria-hidden":"true"},E={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},y={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.307ex",height:"1.91ex",role:"img",focusable:"false",viewBox:"0 -694 1019.6 844","aria-hidden":"true"},u={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},c={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.307ex",height:"1.91ex",role:"img",focusable:"false",viewBox:"0 -694 1019.6 844","aria-hidden":"true"},F={class:"MathJax",jax:"SVG",display:"true",style:{direction:"ltr",display:"block","text-align":"center",margin:"1em 0",position:"relative"}},x={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-1.552ex"},xmlns:"http://www.w3.org/2000/svg",width:"37.823ex",height:"6.801ex",role:"img",focusable:"false",viewBox:"0 -2320 16717.9 3006","aria-hidden":"true"},w={class:"MathJax",jax:"SVG",display:"true",style:{direction:"ltr",display:"block","text-align":"center",margin:"1em 0",position:"relative"}},C={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-1.552ex"},xmlns:"http://www.w3.org/2000/svg",width:"37.823ex",height:"6.801ex",role:"img",focusable:"false",viewBox:"0 -2320 16717.9 3006","aria-hidden":"true"},H={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},f={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.906ex",height:"1.314ex",role:"img",focusable:"false",viewBox:"0 -431 842.6 581","aria-hidden":"true"},V={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},M={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.906ex",height:"1.314ex",role:"img",focusable:"false",viewBox:"0 -431 842.6 581","aria-hidden":"true"},b={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},L={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"8.188ex",height:"1.846ex",role:"img",focusable:"false",viewBox:"0 -666 3619.1 816","aria-hidden":"true"},v={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},B={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"8.188ex",height:"1.846ex",role:"img",focusable:"false",viewBox:"0 -666 3619.1 816","aria-hidden":"true"},D={tabindex:"0"},Z={style:{"text-align":"right"}},A={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},_={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.307ex",height:"1.91ex",role:"img",focusable:"false",viewBox:"0 -694 1019.6 844","aria-hidden":"true"},j={style:{"text-align":"right"}},S={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},P={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.307ex",height:"1.91ex",role:"img",focusable:"false",viewBox:"0 -694 1019.6 844","aria-hidden":"true"},R={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},I={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.049ex",height:"1.934ex",role:"img",focusable:"false",viewBox:"0 -705 905.6 855","aria-hidden":"true"},q={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},N={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.339ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.049ex",height:"1.934ex",role:"img",focusable:"false",viewBox:"0 -705 905.6 855","aria-hidden":"true"},O={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},U={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-2.425ex"},xmlns:"http://www.w3.org/2000/svg",width:"9.211ex",height:"5.981ex",role:"img",focusable:"false",viewBox:"0 -1571.9 4071.1 2643.7","aria-hidden":"true"},G={class:"MathJax",jax:"SVG",display:"true",style:{direction:"ltr",display:"block","text-align":"center",margin:"1em 0",position:"relative"}},z={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-3.275ex"},xmlns:"http://www.w3.org/2000/svg",width:"46.572ex",height:"7.681ex",role:"img",focusable:"false",viewBox:"0 -1947.5 20585 3395","aria-hidden":"true"};function J(X,a,W,Y,$,K){return i(),t("div",null,[a[50]||(a[50]=s("h1",{id:"Bayesian-Updating",tabindex:"-1"},[l("Bayesian Updating "),s("a",{class:"header-anchor",href:"#Bayesian-Updating","aria-label":'Permalink to "Bayesian Updating {#Bayesian-Updating}"'},"​")],-1)),a[51]||(a[51]=s("h2",{id:"Inverse-eigenvalue-problem",tabindex:"-1"},[l("Inverse eigenvalue problem "),s("a",{class:"header-anchor",href:"#Inverse-eigenvalue-problem","aria-label":'Permalink to "Inverse eigenvalue problem {#Inverse-eigenvalue-problem}"'},"​")],-1)),a[52]||(a[52]=s("p",null,"The inverse eigenvalue problem is a classic engineering example. Here we will use Bayesian updating to sample from a bivariate posterior distribution describing unknown quantities of a matrix",-1)),s("mjx-container",r,[(i(),t("svg",Q,a[0]||(a[0]=[n("",1)]))),a[1]||(a[1]=s("mjx-assistive-mml",{unselectable:"on",display:"block",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",overflow:"hidden",width:"100%"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},[s("mrow",{"data-mjx-texclass":"INNER"},[s("mo",{"data-mjx-texclass":"OPEN"},"["),s("mtable",{columnspacing:"1em",rowspacing:"4pt"},[s("mtr",null,[s("mtd",null,[s("msub",null,[s("mi",null,"θ"),s("mn",null,"1")]),s("mo",null,"+"),s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")])]),s("mtd",null,[s("mo",null,"−"),s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")])])]),s("mtr",null,[s("mtd",null,[s("mo",null,"−"),s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")])]),s("mtd",null,[s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")])])])]),s("mo",{"data-mjx-texclass":"CLOSE"},"]")])])],-1))]),s("p",null,[a[6]||(a[6]=l("A matrix of this form can represent different problems, like the stiffness matrix describing a tuned mass damper system. In this example we assume the fixed values ")),s("mjx-container",k,[(i(),t("svg",o,a[2]||(a[2]=[n("",1)]))),a[3]||(a[3]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"θ"),s("mn",null,"1")]),s("mo",null,"="),s("mn",null,"0.5")])],-1))]),a[7]||(a[7]=l(" and ")),s("mjx-container",m,[(i(),t("svg",g,a[4]||(a[4]=[n("",1)]))),a[5]||(a[5]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")]),s("mo",null,"="),s("mn",null,"1.5")])],-1))]),a[8]||(a[8]=l(" for the variables."))]),s("p",null,[a[13]||(a[13]=l("The eigenvalues ")),s("mjx-container",E,[(i(),t("svg",y,a[9]||(a[9]=[n("",1)]))),a[10]||(a[10]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"λ"),s("mn",null,"1")])])],-1))]),a[14]||(a[14]=l(" and ")),s("mjx-container",u,[(i(),t("svg",c,a[11]||(a[11]=[n("",1)]))),a[12]||(a[12]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"λ"),s("mn",null,"2")])])],-1))]),a[15]||(a[15]=l(' of this matrix represent a physical measurable property corrupted by "noise" created for example due to environmental factors or measurement inaccuracy.'))]),s("mjx-container",F,[(i(),t("svg",x,a[16]||(a[16]=[n("",1)]))),a[17]||(a[17]=s("mjx-assistive-mml",{unselectable:"on",display:"block",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",overflow:"hidden",width:"100%"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},[s("msubsup",null,[s("mi",null,"λ"),s("mn",null,"1"),s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",null,"n"),s("mi",null,"o"),s("mi",null,"i"),s("mi",null,"s"),s("mi",null,"y")])]),s("mo",null,"="),s("mfrac",null,[s("mrow",null,[s("mo",{stretchy:"false"},"("),s("msub",null,[s("mi",null,"θ"),s("mn",null,"1")]),s("mo",null,"+"),s("mn",null,"2"),s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")]),s("mo",{stretchy:"false"},")"),s("mo",null,"+"),s("msqrt",null,[s("msubsup",null,[s("mi",null,"θ"),s("mn",null,"1"),s("mn",null,"2")]),s("mo",null,"+"),s("mn",null,"4"),s("msup",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")])]),s("mn",null,"2")])])]),s("mn",null,"2")]),s("mo",null,"+"),s("msub",null,[s("mi",null,"ϵ"),s("mn",null,"1")])])],-1))]),s("mjx-container",w,[(i(),t("svg",C,a[18]||(a[18]=[n("",1)]))),a[19]||(a[19]=s("mjx-assistive-mml",{unselectable:"on",display:"block",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",overflow:"hidden",width:"100%"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},[s("msubsup",null,[s("mi",null,"λ"),s("mn",null,"2"),s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",null,"n"),s("mi",null,"o"),s("mi",null,"i"),s("mi",null,"s"),s("mi",null,"y")])]),s("mo",null,"="),s("mfrac",null,[s("mrow",null,[s("mo",{stretchy:"false"},"("),s("msub",null,[s("mi",null,"θ"),s("mn",null,"1")]),s("mo",null,"+"),s("mn",null,"2"),s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")]),s("mo",{stretchy:"false"},")"),s("mo",null,"−"),s("msqrt",null,[s("msubsup",null,[s("mi",null,"θ"),s("mn",null,"1"),s("mn",null,"2")]),s("mo",null,"+"),s("mn",null,"4"),s("msup",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")])]),s("mn",null,"2")])])]),s("mn",null,"2")]),s("mo",null,"+"),s("msub",null,[s("mi",null,"ϵ"),s("mn",null,"2")])])],-1))]),s("p",null,[a[28]||(a[28]=l('The "noise" terms ')),s("mjx-container",H,[(i(),t("svg",f,a[20]||(a[20]=[n("",1)]))),a[21]||(a[21]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"ϵ"),s("mn",null,"1")])])],-1))]),a[29]||(a[29]=l(" and ")),s("mjx-container",V,[(i(),t("svg",M,a[22]||(a[22]=[n("",1)]))),a[23]||(a[23]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"ϵ"),s("mn",null,"2")])])],-1))]),a[30]||(a[30]=l(" follow a Normal distribution with zero mean and standard deviations ")),s("mjx-container",b,[(i(),t("svg",L,a[24]||(a[24]=[n("",1)]))),a[25]||(a[25]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"σ"),s("mn",null,"1")]),s("mo",null,"="),s("mn",null,"1.0")])],-1))]),a[31]||(a[31]=l(" and ")),s("mjx-container",v,[(i(),t("svg",B,a[26]||(a[26]=[n("",1)]))),a[27]||(a[27]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"σ"),s("mn",null,"2")]),s("mo",null,"="),s("mn",null,"0.1")])],-1))]),a[32]||(a[32]=l("."))]),a[53]||(a[53]=s("p",null,'The synthetic "noisy" data used for the Bayesian updating procedure is given in the following table.',-1)),s("table",D,[s("thead",null,[s("tr",null,[s("th",Z,[s("mjx-container",A,[(i(),t("svg",_,a[33]||(a[33]=[n("",1)]))),a[34]||(a[34]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"λ"),s("mn",null,"1")])])],-1))])]),s("th",j,[s("mjx-container",S,[(i(),t("svg",P,a[35]||(a[35]=[n("",1)]))),a[36]||(a[36]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"λ"),s("mn",null,"2")])])],-1))])])])]),a[37]||(a[37]=s("tbody",null,[s("tr",null,[s("td",{style:{"text-align":"right"}},"1.51"),s("td",{style:{"text-align":"right"}},"0.33")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"4.01"),s("td",{style:{"text-align":"right"}},"0.30")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"3.16"),s("td",{style:{"text-align":"right"}},"0.17")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"3.21"),s("td",{style:{"text-align":"right"}},"0.18")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"2.19"),s("td",{style:{"text-align":"right"}},"0.32")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"1.71"),s("td",{style:{"text-align":"right"}},"0.23")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"2.73"),s("td",{style:{"text-align":"right"}},"0.21")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"5.51"),s("td",{style:{"text-align":"right"}},"0.20")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"1.95"),s("td",{style:{"text-align":"right"}},"0.11")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"4.48"),s("td",{style:{"text-align":"right"}},"0.20")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"1.43"),s("td",{style:{"text-align":"right"}},"0.16")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"2.91"),s("td",{style:{"text-align":"right"}},"0.26")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"3.91"),s("td",{style:{"text-align":"right"}},"0.23")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"3.58"),s("td",{style:{"text-align":"right"}},"0.25")]),s("tr",null,[s("td",{style:{"text-align":"right"}},"2.62"),s("td",{style:{"text-align":"right"}},"0.25")])],-1))]),s("p",null,[a[44]||(a[44]=l("The a priori knowledge of ")),s("mjx-container",R,[(i(),t("svg",I,a[38]||(a[38]=[n("",1)]))),a[39]||(a[39]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"θ"),s("mn",null,"1")])])],-1))]),a[45]||(a[45]=l(" and ")),s("mjx-container",q,[(i(),t("svg",N,a[40]||(a[40]=[n("",1)]))),a[41]||(a[41]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("msub",null,[s("mi",null,"θ"),s("mn",null,"2")])])],-1))]),a[46]||(a[46]=l(" is that they take values between 0.01 and 4. The likelihood function used for this problem is a bivariate Gaussian function with a covariance matrix ")),s("mjx-container",O,[(i(),t("svg",U,a[42]||(a[42]=[n("",1)]))),a[43]||(a[43]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mrow",{"data-mjx-texclass":"INNER"},[s("mo",{"data-mjx-texclass":"OPEN"},"["),s("mtable",{columnspacing:"1em",rowspacing:"4pt"},[s("mtr",null,[s("mtd",null,[s("msubsup",null,[s("mi",null,"σ"),s("mn",null,"1"),s("mn",null,"2")])]),s("mtd",null,[s("mn",null,"0")])]),s("mtr",null,[s("mtd",null,[s("mn",null,"0")]),s("mtd",null,[s("msubsup",null,[s("mi",null,"σ"),s("mn",null,"2"),s("mn",null,"2")])])])]),s("mo",{"data-mjx-texclass":"CLOSE"},"]")])])],-1))]),a[47]||(a[47]=l(", with off-diagonal terms equal to 0 and the diagonal terms corresponding to the variances of the respective noise terms."))]),s("mjx-container",G,[(i(),t("svg",z,a[48]||(a[48]=[n("",1)]))),a[49]||(a[49]=s("mjx-assistive-mml",{unselectable:"on",display:"block",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",overflow:"hidden",width:"100%"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},[s("mi",null,"P"),s("mo",{stretchy:"false"},"("),s("mi",null,"λ"),s("mo",{"data-mjx-texclass":"ORD",stretchy:"false"},"|"),s("mi",null,"θ"),s("mo",{stretchy:"false"},")"),s("mo",null,"∝"),s("mi",null,"exp"),s("mo",{"data-mjx-texclass":"NONE"},"⁡"),s("mrow",{"data-mjx-texclass":"INNER"},[s("mo",{"data-mjx-texclass":"OPEN"},"["),s("mo",null,"−"),s("mfrac",null,[s("mn",null,"1"),s("mn",null,"2")]),s("munderover",null,[s("mo",{"data-mjx-texclass":"OP"},"∑"),s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",null,"i"),s("mo",null,"="),s("mn",null,"1")]),s("mn",null,"2")]),s("munderover",null,[s("mo",{"data-mjx-texclass":"OP"},"∑"),s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",null,"n"),s("mo",null,"="),s("mn",null,"1")]),s("mrow",{"data-mjx-texclass":"ORD"},[s("mn",null,"15")])]),s("msup",null,[s("mrow",{"data-mjx-texclass":"ORD"},[s("mrow",{"data-mjx-texclass":"INNER"},[s("mo",{"data-mjx-texclass":"OPEN"},"("),s("mfrac",null,[s("mrow",null,[s("msubsup",null,[s("mi",null,"λ"),s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",null,"i"),s("mo",null,","),s("mi",null,"n")]),s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",null,"d"),s("mi",null,"a"),s("mi",null,"t"),s("mi",null,"a")])]),s("mo",null,"−"),s("msubsup",null,[s("mi",null,"λ"),s("mi",null,"i"),s("mrow",{"data-mjx-texclass":"ORD"},[s("mi",null,"m"),s("mi",null,"o"),s("mi",null,"d"),s("mi",null,"e"),s("mi",null,"l")])])]),s("msub",null,[s("mi",null,"σ"),s("mi",null,"i")])]),s("mo",{"data-mjx-texclass":"CLOSE"},")")])]),s("mn",null,"2")]),s("mo",{"data-mjx-texclass":"CLOSE"},"]")])])],-1))]),a[54]||(a[54]=n("",27))])}const t1=e(d,[["render",J]]);export{a1 as __pageData,t1 as default};
