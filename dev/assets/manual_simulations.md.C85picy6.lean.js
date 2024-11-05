import{_ as t,c as a,a5 as i,o as s}from"./chunks/framework.Dt0Z5R5R.js";const n="/UncertaintyQuantification.jl/dev/assets/faure-sequence.D7aVGsOn.svg",m=JSON.parse('{"title":"Simulations","description":"","frontmatter":{},"headers":[],"relativePath":"manual/simulations.md","filePath":"manual/simulations.md","lastUpdated":null}'),l={name:"manual/simulations.md"};function o(d,e,h,r,p,c){return s(),a("div",null,e[0]||(e[0]=[i(`<h1 id="simulations" tabindex="-1">Simulations <a class="header-anchor" href="#simulations" aria-label="Permalink to &quot;Simulations&quot;">​</a></h1><h2 id="Monte-Carlo" tabindex="-1">Monte Carlo <a class="header-anchor" href="#Monte-Carlo" aria-label="Permalink to &quot;Monte Carlo {#Monte-Carlo}&quot;">​</a></h2><p>The Monte-Carlo (MC) method is a method of sampling random numbers that dates back to 1777. The name was suggested by Nick Metropolis when MC was used while working on the Manhattan Project. It is used in Random Number Generators which generally produce pseudo random numbers.</p><h2 id="Quasi-Monte-Carlo" tabindex="-1">Quasi Monte Carlo <a class="header-anchor" href="#Quasi-Monte-Carlo" aria-label="Permalink to &quot;Quasi Monte Carlo {#Quasi-Monte-Carlo}&quot;">​</a></h2><p>Quasi Monte Carlo (QMC), is a method of producing samples similar to those generated via Monte Carlo (MC). The difference being that QMC samples are generated deterministically in a way to ensure they are evenly distributed across the sampling space, not forming clutters or voids as MC samples might. This makes QMC more efficient than MC for lots of applications since fewer samples are needed in order to produce a sufficient density of samples throughout. There are multiple ways of QMC-sampling which can be classified as either digital nets or lattices. [<a href="/UncertaintyQuantification.jl/dev/references#owenQuasiMonteCarlo2009">16</a>]</p><p>Included here are <code>LatticeRuleSampling</code> and the digital nets <code>SobolSampling</code>, <code>HaltonSampling</code>, <code>FaureSampling</code> and <code>LatinHaypercubeSampling</code>.</p><p>However, being deterministic, these QMC samples are missing the properties related to randomness that MC samples have. To gain these properties it is possible to randomize QMC samples. There are several randomization methods, useful in different cases, depending on the QMC method in use. [<a href="/UncertaintyQuantification.jl/dev/references#owenQuasiMonteCarlo2009">16</a>]</p><p>Implemented in this package are Owen-Scramble and Matousek-Scramble, two similar methods useful for Sobol and Faure Sampling aswell as Shift which can be used for Lattice Rule Sampling. There also is an algorithm for Halton Sampling, that constructs builds samples from the ground up as opposed to randomizing existing samples which is what the aforementioned methods do. [<a href="/UncertaintyQuantification.jl/dev/references#owenRandomizedHalton2017">17</a>]</p><p>To sample using one of these methods, simply create an instance of the corresponding struct with the desired parameters and then call the sample function with this instance. The parameters are <code>n::Integer</code> which is the number of samples, and <code>randomization::Symbol</code> which encodes the randomization method that should be used. The different possible symbols are: <code>:none</code>, <code>:matousek</code>, <code>:owen</code>, <code>:shift</code> and <code>:randomizedhalton</code>.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> RandomVariable</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Uniform</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    qmc </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> LatinHypercubeSampling</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    samples </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sample</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, qmc)</span></span></code></pre></div><p>Note that not all randomization methods are possible to use for every QMC-method. Also, if no <code>randomization</code>-symbol is given, the default will be used. View the following table for details.</p><table tabindex="0"><thead><tr><th style="text-align:left;">QMC-method</th><th style="text-align:left;">DEFAULT</th><th style="text-align:center;">:matousek</th><th style="text-align:center;">:owen</th><th style="text-align:center;">:shift</th><th style="text-align:center;">:randomizedhalton</th><th style="text-align:center;">:none</th></tr></thead><tbody><tr><td style="text-align:left;">LatticeRuleSampling</td><td style="text-align:left;">:shift</td><td style="text-align:center;">❌</td><td style="text-align:center;">❌</td><td style="text-align:center;">✅</td><td style="text-align:center;">❌</td><td style="text-align:center;">✅</td></tr><tr><td style="text-align:left;">SobolSampling</td><td style="text-align:left;">:matousek</td><td style="text-align:center;">✅</td><td style="text-align:center;">✅</td><td style="text-align:center;">❌</td><td style="text-align:center;">❌</td><td style="text-align:center;">✅</td></tr><tr><td style="text-align:left;">FaureSampling</td><td style="text-align:left;">:matousek</td><td style="text-align:center;">✅</td><td style="text-align:center;">✅</td><td style="text-align:center;">❌</td><td style="text-align:center;">❌</td><td style="text-align:center;">✅</td></tr><tr><td style="text-align:left;">HaltonSampling</td><td style="text-align:left;">:randomizedhalton</td><td style="text-align:center;">❌</td><td style="text-align:center;">❌</td><td style="text-align:center;">❌</td><td style="text-align:center;">✅</td><td style="text-align:center;">✅</td></tr></tbody></table><div class="tip custom-block"><p class="custom-block-title">Note</p><p><code>LatinHypercubeSampling</code> is already random and thus doesn&#39;t have the <code>randomization</code> parameter.</p></div><p>It is of course possible to directly create the struct inside the <code>sample</code>-call, enabling a more efficient version of the example above which looks like this:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> RandomVariable</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Uniform</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    samples </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> sample</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">LatinHypercubeSampling</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><div class="tip custom-block"><p class="custom-block-title">Note</p><p>When chosing <code>n</code>, bear in mind that for <code>SobolSampling</code> and <code>FaureSampling</code>, <code>n</code> must fit the base that is used for creating the respective sequence. For <code>SobolSampling</code> the base is always equal to 2 while for <code>FaureSampling</code>, it depends on the number of input-variables. If <code>n</code> is not a power of the base, it will automatically be increased to the next power.</p></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> RandomVariable</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Uniform</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    samples </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> SobolSampling</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>SobolSampling(128, :matousek)</span></span></code></pre></div><p>To emphasize the importance of randomization, look at the correlations that might occur using unrandomized qmc and how they are fixed by randomizing.</p><p>This is the 7th dimension plotted against the 8th in Faure Sampling, unrandomized vs. randomized via Owen Scramble:</p><p><img src="`+n+'" alt=""></p>',21)]))}const g=t(l,[["render",o]]);export{m as __pageData,g as default};
