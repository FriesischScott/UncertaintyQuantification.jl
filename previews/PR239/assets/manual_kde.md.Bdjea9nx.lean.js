import{_ as l,c as Q,o as s,j as t,ai as i,a}from"./chunks/framework.XrgkmTJE.js";const n="/UncertaintyQuantification.jl/previews/PR239/assets/kernel-density.BWQbJZ-p.svg",M=JSON.parse('{"title":"Kernel Density Estimation","description":"","frontmatter":{},"headers":[],"relativePath":"manual/kde.md","filePath":"manual/kde.md","lastUpdated":null}'),T={name:"manual/kde.md"},o={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},r={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.464ex"},xmlns:"http://www.w3.org/2000/svg",width:"2.218ex",height:"2.891ex",role:"img",focusable:"false",viewBox:"0 -1073 980.3 1278","aria-hidden":"true"},d={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},m={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.439ex"},xmlns:"http://www.w3.org/2000/svg",width:"10.923ex",height:"1.984ex",role:"img",focusable:"false",viewBox:"0 -683 4827.8 877","aria-hidden":"true"},h={class:"MathJax",jax:"SVG",display:"true",style:{direction:"ltr",display:"block","text-align":"center",margin:"1em 0",position:"relative"}},p={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-2.819ex"},xmlns:"http://www.w3.org/2000/svg",width:"33.515ex",height:"6.354ex",role:"img",focusable:"false",viewBox:"0 -1562.5 14813.6 2808.5","aria-hidden":"true"},g={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},u={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.025ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.303ex",height:"1.595ex",role:"img",focusable:"false",viewBox:"0 -694 576 705","aria-hidden":"true"},k={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},x={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"0"},xmlns:"http://www.w3.org/2000/svg",width:"2.011ex",height:"1.545ex",role:"img",focusable:"false",viewBox:"0 -683 889 683","aria-hidden":"true"},c={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},w={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.025ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.303ex",height:"1.595ex",role:"img",focusable:"false",viewBox:"0 -694 576 705","aria-hidden":"true"};function y(f,e,L,H,b,E){return s(),Q("div",null,[e[29]||(e[29]=t("h1",{id:"Kernel-Density-Estimation",tabindex:"-1"},[a("Kernel Density Estimation "),t("a",{class:"header-anchor",href:"#Kernel-Density-Estimation","aria-label":'Permalink to "Kernel Density Estimation {#Kernel-Density-Estimation}"'},"​")],-1)),e[30]||(e[30]=t("p",null,[a("Kernel density estimation (KDE) is a non-parametric method to estimate the probability density function of a random variable through "),t("em",null,"kernel smoothing"),a(" ["),t("a",{href:"/UncertaintyQuantification.jl/previews/PR239/references#silvermanDensityEstimationStatistics1986"},"9"),a("].")],-1)),t("p",null,[e[4]||(e[4]=a("The kernel density estimate ")),t("mjx-container",o,[(s(),Q("svg",r,e[0]||(e[0]=[i("",1)]))),e[1]||(e[1]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("msub",null,[t("mrow",{"data-mjx-texclass":"ORD"},[t("mover",null,[t("mi",null,"f"),t("mo",{stretchy:"false"},"^")])]),t("mi",null,"h")])])],-1))]),e[5]||(e[5]=a(" of a univariate density ")),e[6]||(e[6]=t("code",null,"f",-1)),e[7]||(e[7]=a(" based on a random sample ")),t("mjx-container",d,[(s(),Q("svg",m,e[2]||(e[2]=[i("",1)]))),e[3]||(e[3]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("msub",null,[t("mi",null,"X"),t("mn",null,"1")]),t("mo",null,","),t("mo",null,"…"),t("mo",null,","),t("msub",null,[t("mi",null,"X"),t("mi",null,"n")])])],-1))]),e[8]||(e[8]=a(" is defined as"))]),t("mjx-container",h,[(s(),Q("svg",p,e[9]||(e[9]=[i("",1)]))),e[10]||(e[10]=t("mjx-assistive-mml",{unselectable:"on",display:"block",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",overflow:"hidden",width:"100%"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},[t("msub",null,[t("mrow",{"data-mjx-texclass":"ORD"},[t("mover",null,[t("mi",null,"f"),t("mo",{stretchy:"false"},"^")])]),t("mi",null,"h")]),t("mo",{stretchy:"false"},"("),t("mi",null,"x"),t("mo",{stretchy:"false"},")"),t("mo",null,"="),t("msup",null,[t("mi",null,"n"),t("mrow",{"data-mjx-texclass":"ORD"},[t("mo",null,"−"),t("mn",null,"1")])]),t("munderover",null,[t("mo",{"data-mjx-texclass":"OP"},"∑"),t("mrow",{"data-mjx-texclass":"ORD"},[t("mi",null,"i"),t("mo",null,"="),t("mn",null,"1")]),t("mi",null,"n")]),t("msup",null,[t("mi",null,"h"),t("mrow",{"data-mjx-texclass":"ORD"},[t("mo",null,"−"),t("mn",null,"1")])]),t("mi",null,"K"),t("mrow",{"data-mjx-texclass":"INNER"},[t("mo",{"data-mjx-texclass":"OPEN"},"{"),t("mfrac",null,[t("mrow",null,[t("mi",null,"x"),t("mo",null,"−"),t("msub",null,[t("mi",null,"X"),t("mi",null,"i")])]),t("mi",null,"h")]),t("mo",{"data-mjx-texclass":"CLOSE"},"}")]),t("mo",null,",")])],-1))]),t("p",null,[e[17]||(e[17]=a("where ")),t("mjx-container",g,[(s(),Q("svg",u,e[11]||(e[11]=[t("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[t("g",{"data-mml-node":"math"},[t("g",{"data-mml-node":"mi"},[t("path",{"data-c":"210E",d:"M137 683Q138 683 209 688T282 694Q294 694 294 685Q294 674 258 534Q220 386 220 383Q220 381 227 388Q288 442 357 442Q411 442 444 415T478 336Q478 285 440 178T402 50Q403 36 407 31T422 26Q450 26 474 56T513 138Q516 149 519 151T535 153Q555 153 555 145Q555 144 551 130Q535 71 500 33Q466 -10 419 -10H414Q367 -10 346 17T325 74Q325 90 361 192T398 345Q398 404 354 404H349Q266 404 205 306L198 293L164 158Q132 28 127 16Q114 -11 83 -11Q69 -11 59 -2T48 16Q48 30 121 320L195 616Q195 629 188 632T149 637H128Q122 643 122 645T124 664Q129 683 137 683Z",style:{"stroke-width":"3"}})])])],-1)]))),e[12]||(e[12]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("mi",null,"h")])],-1))]),e[18]||(e[18]=a(" is the so called ")),e[19]||(e[19]=t("em",null,"bandwidth",-1)),e[20]||(e[20]=a(" and ")),t("mjx-container",k,[(s(),Q("svg",x,e[13]||(e[13]=[t("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[t("g",{"data-mml-node":"math"},[t("g",{"data-mml-node":"mi"},[t("path",{"data-c":"1D43E",d:"M285 628Q285 635 228 637Q205 637 198 638T191 647Q191 649 193 661Q199 681 203 682Q205 683 214 683H219Q260 681 355 681Q389 681 418 681T463 682T483 682Q500 682 500 674Q500 669 497 660Q496 658 496 654T495 648T493 644T490 641T486 639T479 638T470 637T456 637Q416 636 405 634T387 623L306 305Q307 305 490 449T678 597Q692 611 692 620Q692 635 667 637Q651 637 651 648Q651 650 654 662T659 677Q662 682 676 682Q680 682 711 681T791 680Q814 680 839 681T869 682Q889 682 889 672Q889 650 881 642Q878 637 862 637Q787 632 726 586Q710 576 656 534T556 455L509 418L518 396Q527 374 546 329T581 244Q656 67 661 61Q663 59 666 57Q680 47 717 46H738Q744 38 744 37T741 19Q737 6 731 0H720Q680 3 625 3Q503 3 488 0H478Q472 6 472 9T474 27Q478 40 480 43T491 46H494Q544 46 544 71Q544 75 517 141T485 216L427 354L359 301L291 248L268 155Q245 63 245 58Q245 51 253 49T303 46H334Q340 37 340 35Q340 19 333 5Q328 0 317 0Q314 0 280 1T180 2Q118 2 85 2T49 1Q31 1 31 11Q31 13 34 25Q38 41 42 43T65 46Q92 46 125 49Q139 52 144 61Q147 65 216 339T285 628Z",style:{"stroke-width":"3"}})])])],-1)]))),e[14]||(e[14]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("mi",null,"K")])],-1))]),e[21]||(e[21]=a(" is the kernel function. The kernel function is assumed to be a symmetric probability density and is set to be a Gaussian density in ")),e[22]||(e[22]=t("em",null,"UncertaintyQuantification.jl",-1)),e[23]||(e[23]=a(". The bandwidth ")),t("mjx-container",c,[(s(),Q("svg",w,e[15]||(e[15]=[t("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[t("g",{"data-mml-node":"math"},[t("g",{"data-mml-node":"mi"},[t("path",{"data-c":"210E",d:"M137 683Q138 683 209 688T282 694Q294 694 294 685Q294 674 258 534Q220 386 220 383Q220 381 227 388Q288 442 357 442Q411 442 444 415T478 336Q478 285 440 178T402 50Q403 36 407 31T422 26Q450 26 474 56T513 138Q516 149 519 151T535 153Q555 153 555 145Q555 144 551 130Q535 71 500 33Q466 -10 419 -10H414Q367 -10 346 17T325 74Q325 90 361 192T398 345Q398 404 354 404H349Q266 404 205 306L198 293L164 158Q132 28 127 16Q114 -11 83 -11Q69 -11 59 -2T48 16Q48 30 121 320L195 616Q195 629 188 632T149 637H128Q122 643 122 645T124 664Q129 683 137 683Z",style:{"stroke-width":"3"}})])])],-1)]))),e[16]||(e[16]=t("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[t("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[t("mi",null,"h")])],-1))]),e[24]||(e[24]=a(" also called the ")),e[25]||(e[25]=t("em",null,"smoothing parameter",-1)),e[26]||(e[26]=a(" has a strong effect on the resulting density estimate. There are various different methods to select an optimal bandwith. Here we have decided to apply the method developed by Sheather & Jones [")),e[27]||(e[27]=t("a",{href:"/UncertaintyQuantification.jl/previews/PR239/references#sheatherReliableDataBasedBandwidth1991"},"10",-1)),e[28]||(e[28]=a("] for its excellent performance and straightforward implementation."))]),e[31]||(e[31]=i("",8))])}const D=l(T,[["render",y]]);export{M as __pageData,D as default};
