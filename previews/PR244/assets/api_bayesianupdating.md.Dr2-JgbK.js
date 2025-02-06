import{_ as o,C as l,c as p,o as r,ai as t,j as a,a as s,G as n}from"./chunks/framework.B6ywfKPV.js";const U=JSON.parse('{"title":"Bayesian Updating","description":"","frontmatter":{},"headers":[],"relativePath":"api/bayesianupdating.md","filePath":"api/bayesianupdating.md","lastUpdated":null}'),d={name:"api/bayesianupdating.md"},h={class:"jldocstring custom-block",open:""},c={class:"jldocstring custom-block",open:""},u={class:"jldocstring custom-block",open:""},k={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""},y={class:"jldocstring custom-block",open:""},b={class:"jldocstring custom-block",open:""};function m(f,i,E,C,v,F){const e=l("Badge");return r(),p("div",null,[i[21]||(i[21]=t('<h1 id="Bayesian-Updating" tabindex="-1">Bayesian Updating <a class="header-anchor" href="#Bayesian-Updating" aria-label="Permalink to &quot;Bayesian Updating {#Bayesian-Updating}&quot;">​</a></h1><p>Methods for Bayesian updating.</p><h2 id="index" tabindex="-1">Index <a class="header-anchor" href="#index" aria-label="Permalink to &quot;Index&quot;">​</a></h2><ul><li><a href="#UncertaintyQuantification.AbstractBayesianMethod"><code>UncertaintyQuantification.AbstractBayesianMethod</code></a></li><li><a href="#UncertaintyQuantification.AbstractBayesianPointEstimate"><code>UncertaintyQuantification.AbstractBayesianPointEstimate</code></a></li><li><a href="#UncertaintyQuantification.MaximumAPosterioriBayesian"><code>UncertaintyQuantification.MaximumAPosterioriBayesian</code></a></li><li><a href="#UncertaintyQuantification.MaximumLikelihoodBayesian"><code>UncertaintyQuantification.MaximumLikelihoodBayesian</code></a></li><li><a href="#UncertaintyQuantification.SingleComponentMetropolisHastings"><code>UncertaintyQuantification.SingleComponentMetropolisHastings</code></a></li><li><a href="#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo"><code>UncertaintyQuantification.TransitionalMarkovChainMonteCarlo</code></a></li><li><a href="#UncertaintyQuantification.bayesianupdating"><code>UncertaintyQuantification.bayesianupdating</code></a></li></ul><h2 id="types" tabindex="-1">Types <a class="header-anchor" href="#types" aria-label="Permalink to &quot;Types&quot;">​</a></h2>',5)),a("details",h,[a("summary",null,[i[0]||(i[0]=a("a",{id:"UncertaintyQuantification.AbstractBayesianMethod",href:"#UncertaintyQuantification.AbstractBayesianMethod"},[a("span",{class:"jlbinding"},"UncertaintyQuantification.AbstractBayesianMethod")],-1)),i[1]||(i[1]=s()),n(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[2]||(i[2]=t('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">AbstractBayesianMethod</span></span></code></pre></div><p>Subtypes are used to dispatch to the differenct MCMC methods in <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a>.</p><p>Subtypes are:</p><ul><li><p><a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.SingleComponentMetropolisHastings"><code>SingleComponentMetropolisHastings</code></a></p></li><li><p><a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo"><code>TransitionalMarkovChainMonteCarlo</code></a></p></li></ul><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/57c46daf5a8cbb2d918b6b53ff9d3f2ca9d09851/src/UncertaintyQuantification.jl#L46-L55" target="_blank" rel="noreferrer">source</a></p>',5))]),a("details",c,[a("summary",null,[i[3]||(i[3]=a("a",{id:"UncertaintyQuantification.SingleComponentMetropolisHastings",href:"#UncertaintyQuantification.SingleComponentMetropolisHastings"},[a("span",{class:"jlbinding"},"UncertaintyQuantification.SingleComponentMetropolisHastings")],-1)),i[4]||(i[4]=s()),n(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[5]||(i[5]=t('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SingleComponentMetropolisHastings</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(proposal, x0, n, burnin, islog)</span></span></code></pre></div><p>Passed to <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a> to run the single-component Metropolis-Hastings algorithm starting from <code>x0</code> with univariate proposal distibution <code>proposal</code>. Will generate <code>n</code> samples <em>after</em> performing <code>burnin</code> steps of the Markov chain and discarding the samples. The flag <code>islog</code> specifies whether the prior and likelihood functions passed to the <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a> method are already given as logarithms.</p><p>Alternative constructor</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    SingleComponentMetropolisHastings</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(proposal, x0, n, burnin)  </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># `islog` = true</span></span></code></pre></div><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/57c46daf5a8cbb2d918b6b53ff9d3f2ca9d09851/src/modelupdating/bayesianupdating.jl#L1-L12" target="_blank" rel="noreferrer">source</a></p>',5))]),a("details",u,[a("summary",null,[i[6]||(i[6]=a("a",{id:"UncertaintyQuantification.TransitionalMarkovChainMonteCarlo",href:"#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo"},[a("span",{class:"jlbinding"},"UncertaintyQuantification.TransitionalMarkovChainMonteCarlo")],-1)),i[7]||(i[7]=s()),n(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[8]||(i[8]=t('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">TransitionalMarkovChainMonteCarlo</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, n, burnin, β, islog)</span></span>\n<span class="line"></span>\n<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Passed to [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">`bayesianupdating`</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">](</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@ref</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) to run thetransitional Markov chain Monte Carlo algorithm  with [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">`RandomVariable&#39;](@ref) vector `</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">prior</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">`. At each transitional level, one sample will be generated from `</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">n</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">` independent Markov chains after `</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">burnin</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">` steps have been discarded. The flag `</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">islog</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">` specifies whether the prior and likelihood functions passed to the  [`</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">bayesianupdating</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">`](@ref) method are already  given as logarithms.</span></span></code></pre></div><p>Alternative constructors</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    TransitionalMarkovChainMonteCarlo</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, n, burnin, β)  </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># `islog` = true</span></span>\n<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">     TransitionalMarkovChainMonteCarlo</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, n, burnin)    </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># `β` = 0.2,  `islog` = true</span></span></code></pre></div><p><strong>References</strong></p><p>[<a href="/UncertaintyQuantification.jl/previews/PR244/references#chingTransitionalMarkovChain2007">21</a>]</p><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/57c46daf5a8cbb2d918b6b53ff9d3f2ca9d09851/src/modelupdating/bayesianupdating.jl#L134-L150" target="_blank" rel="noreferrer">source</a></p>',6))]),a("details",k,[a("summary",null,[i[9]||(i[9]=a("a",{id:"UncertaintyQuantification.AbstractBayesianPointEstimate",href:"#UncertaintyQuantification.AbstractBayesianPointEstimate"},[a("span",{class:"jlbinding"},"UncertaintyQuantification.AbstractBayesianPointEstimate")],-1)),i[10]||(i[10]=s()),n(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[11]||(i[11]=t('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">AbstractBayesianPointEstimate</span></span></code></pre></div><p>Subtypes are used to dispatch to the differenct point estimation methods in <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a>.</p><p>Subtypes are:</p><ul><li><p><a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.MaximumAPosterioriBayesian"><code>MaximumAPosterioriBayesian</code></a></p></li><li><p><a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.MaximumLikelihoodBayesian"><code>MaximumLikelihoodBayesian</code></a></p></li></ul><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/57c46daf5a8cbb2d918b6b53ff9d3f2ca9d09851/src/UncertaintyQuantification.jl#L58-L67" target="_blank" rel="noreferrer">source</a></p>',5))]),a("details",g,[a("summary",null,[i[12]||(i[12]=a("a",{id:"UncertaintyQuantification.MaximumAPosterioriBayesian",href:"#UncertaintyQuantification.MaximumAPosterioriBayesian"},[a("span",{class:"jlbinding"},"UncertaintyQuantification.MaximumAPosterioriBayesian")],-1)),i[13]||(i[13]=s()),n(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[14]||(i[14]=t('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">MaximumAPosterioriBayesian</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, optimmethod, x0; islog, lowerbounds, upperbounds)</span></span></code></pre></div><p>Passed to <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a> to estimate one or more maxima of the posterior distribution starting from <code>x0</code>. The optimization uses the method specified in <code>optimmethod</code>. Will calculate one estimation per point in x0. The flag <code>islog</code> specifies whether the prior and likelihood functions passed to the <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a> method are already given as logarithms. <code>lowerbounds</code> and <code>upperbounds</code> specify optimization intervals.</p><p>Alternative constructors</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    MaximumAPosterioriBayesian</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, optimmethod, x0; islog) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># `lowerbounds` = [-Inf], # `upperbounds` = [Inf]</span></span>\n<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    MaximumAPosterioriBayesian</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, optimmethod, x0)  </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># `islog` = true</span></span></code></pre></div><p>See also <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.MaximumLikelihoodBayesian"><code>MaximumLikelihoodBayesian</code></a>, <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a>, <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo"><code>TransitionalMarkovChainMonteCarlo</code></a>.</p><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/57c46daf5a8cbb2d918b6b53ff9d3f2ca9d09851/src/modelupdating/bayesianMAP.jl#L1-L13" target="_blank" rel="noreferrer">source</a></p>',6))]),a("details",y,[a("summary",null,[i[15]||(i[15]=a("a",{id:"UncertaintyQuantification.MaximumLikelihoodBayesian",href:"#UncertaintyQuantification.MaximumLikelihoodBayesian"},[a("span",{class:"jlbinding"},"UncertaintyQuantification.MaximumLikelihoodBayesian")],-1)),i[16]||(i[16]=s()),n(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[17]||(i[17]=t('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">MaximumLikelihoodBayesian</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, optimmethod, x0; islog, lowerbounds, upperbounds)</span></span></code></pre></div><p>Passed to <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a> to estimate one or more maxima of the likelihood starting from <code>x0</code>. The optimization uses the method specified in <code>optimmethod</code>. Will calculate one estimation per point in x0. The flag <code>islog</code> specifies whether the prior and likelihood functions passed to the <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a> method are already given as logarithms. <code>lowerbounds</code> and <code>upperbounds</code> specify optimization intervals.</p><p>Alternative constructors</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    MaximumLikelihoodBayesian</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, optimmethod, x0; islog) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># `lowerbounds` = [-Inf], # `upperbounds` = [Inf]</span></span>\n<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    MaximumLikelihoodBayesian</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, optimmethod, x0)  </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># `islog` = true</span></span></code></pre></div><p><strong>Notes</strong></p><p>The method uses <code>prior</code> only as information on which parameters are supposed to be optimized. The prior itself does not influence the result of the maximum likelihood estimate and can be given as a dummy distribution. For example, if two parameters <code>a</code> and <code>b</code> are supposed to be optimized, the prior could look like this</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    prior </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> RandomVariable</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Uniform</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:a</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:b</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span></code></pre></div><p>See also <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.MaximumAPosterioriBayesian"><code>MaximumAPosterioriBayesian</code></a>, <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a>, <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo"><code>TransitionalMarkovChainMonteCarlo</code></a>.</p><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/57c46daf5a8cbb2d918b6b53ff9d3f2ca9d09851/src/modelupdating/bayesianMAP.jl#L46-L64" target="_blank" rel="noreferrer">source</a></p>',9))]),i[22]||(i[22]=a("h2",{id:"methods",tabindex:"-1"},[s("Methods "),a("a",{class:"header-anchor",href:"#methods","aria-label":'Permalink to "Methods"'},"​")],-1)),a("details",b,[a("summary",null,[i[18]||(i[18]=a("a",{id:"UncertaintyQuantification.bayesianupdating",href:"#UncertaintyQuantification.bayesianupdating"},[a("span",{class:"jlbinding"},"UncertaintyQuantification.bayesianupdating")],-1)),i[19]||(i[19]=s()),n(e,{type:"info",class:"jlObjectType jlFunction",text:"Function"})]),i[20]||(i[20]=t(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">bayesianupdating</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(likelihood, models, pointestimate; prior)</span></span></code></pre></div><p>Perform bayesian updating using the given <code>likelihood</code>, <code>models</code> and any point estimation method <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.AbstractBayesianPointEstimate"><code>AbstractBayesianPointEstimate</code></a>.</p><p><strong>Notes</strong></p><p>Method can be called with an empty Vector of models, i.e.</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>bayesianupdating(likelihood, [], pointestimate)</span></span></code></pre></div><p>If <code>prior</code> is not given, the method will construct a prior distribution from the prior specified in <code>AbstractBayesianPointEstimate.prior</code>.</p><p><code>likelihood</code> is a Julia function which must be defined in terms of a <code>DataFrame</code> of samples, and must evaluate the likelihood for each row of the <code>DataFrame</code></p><p>For example, a loglikelihood based on normal distribution using &#39;Data&#39;:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">likelihood</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(df) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sum</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">logpdf</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Normal</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(df_i</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">x, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), Data)) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> df_i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> eachrow</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(df)]</span></span></code></pre></div><p>If a model evaluation is required to evaluate the likelihood, a vector of <code>UQModel</code>s must be passed to <code>bayesianupdating</code>. For example if the variable <code>x</code> above is the output of a numerical model.</p><p>For a general overview of the function, see <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.bayesianupdating"><code>bayesianupdating</code></a>.</p><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/57c46daf5a8cbb2d918b6b53ff9d3f2ca9d09851/src/modelupdating/bayesianMAP.jl#L100-L124" target="_blank" rel="noreferrer">source</a></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">bayesianupdating</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prior, likelihood, models, mcmc)</span></span></code></pre></div><p>Perform bayesian updating using the given <code>prior</code>, <code>likelihood</code>, <code>models</code> and any MCMC sampler <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.AbstractBayesianMethod"><code>AbstractBayesianMethod</code></a>.</p><p>Alternatively the method can be called without <code>models</code>.</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>bayesianupdating(prior, likelihood, mcmc)</span></span></code></pre></div><p>When using <a href="/UncertaintyQuantification.jl/previews/PR244/api/bayesianupdating#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo"><code>TransitionalMarkovChainMonteCarlo</code></a> the <code>prior</code> can automatically be constructed.</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>bayesinupdating(likelihood, models, tmcmc)</span></span>
<span class="line"><span>bayesianupdating(likelihood, tmcmc)</span></span></code></pre></div><p><strong>Notes</strong></p><p><code>likelihood</code> is a Julia function which must be defined in terms of a <code>DataFrame</code> of samples, and must evaluate the likelihood for each row of the <code>DataFrame</code></p><p>For example, a loglikelihood based on normal distribution using &#39;Data&#39;:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">likelihood</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(df) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sum</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">logpdf</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Normal</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(df_i</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">x, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), Data)) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> df_i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> eachrow</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(df)]</span></span></code></pre></div><p>If a model evaluation is required to evaluate the likelihood, a vector of <code>UQModel</code>s must be passed to <code>bayesianupdating</code>. For example if the variable <code>x</code> above is the output of a numerical model.</p><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/57c46daf5a8cbb2d918b6b53ff9d3f2ca9d09851/src/modelupdating/bayesianupdating.jl#L34-L61" target="_blank" rel="noreferrer">source</a></p>`,24))])])}const B=o(d,[["render",m]]);export{U as __pageData,B as default};
