import{_ as l,C as h,c as p,o as k,ai as a,j as s,a as e,G as n}from"./chunks/framework.DCiQS8J1.js";const f=JSON.parse('{"title":"Simulations","description":"","frontmatter":{},"headers":[],"relativePath":"api/simulations.md","filePath":"api/simulations.md","lastUpdated":null}'),r={name:"api/simulations.md"},o={class:"jldocstring custom-block",open:""},d={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""};function E(c,i,y,u,b,F){const t=h("Badge");return k(),p("div",null,[i[9]||(i[9]=a('<h1 id="simulations" tabindex="-1">Simulations <a class="header-anchor" href="#simulations" aria-label="Permalink to &quot;Simulations&quot;">​</a></h1><p>Various Monte Carlo based simulations for a wide range of applications.</p><h2 id="index" tabindex="-1">Index <a class="header-anchor" href="#index" aria-label="Permalink to &quot;Index&quot;">​</a></h2><ul><li><a href="#UncertaintyQuantification.SubSetInfinity"><code>UncertaintyQuantification.SubSetInfinity</code></a></li><li><a href="#UncertaintyQuantification.SubSetInfinityAdaptive"><code>UncertaintyQuantification.SubSetInfinityAdaptive</code></a></li><li><a href="#UncertaintyQuantification.SubSetSimulation"><code>UncertaintyQuantification.SubSetSimulation</code></a></li></ul><h2 id="types" tabindex="-1">Types <a class="header-anchor" href="#types" aria-label="Permalink to &quot;Types&quot;">​</a></h2>',5)),s("details",o,[s("summary",null,[i[0]||(i[0]=s("a",{id:"UncertaintyQuantification.SubSetSimulation",href:"#UncertaintyQuantification.SubSetSimulation"},[s("span",{class:"jlbinding"},"UncertaintyQuantification.SubSetSimulation")],-1)),i[1]||(i[1]=e()),n(t,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[2]||(i[2]=a(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SubSetSimulation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(n</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, target</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Float64</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, levels</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, proposal</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">UnivariateDistribution</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Defines the properties of a Subset simulation where <code>n</code> is the number of initial samples, <code>target</code> is the target probability of failure at each level, <code>levels</code> is the maximum number of levels and <code>proposal</code> is the proposal distribution for the markov chain monte carlo.</p><p><strong>Examples</strong></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> SubSetSimulation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Uniform</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SubSetSimulation</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Uniform{Float64}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><p><strong>References</strong></p><p>[<a href="/UncertaintyQuantification.jl/previews/PR243/references#auEstimationSmallFailure2001">14</a>]</p><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/9573d7c5373e6153b9cbb7558b521a09977dfdb7/src/simulations/subset.jl#L3-L20" target="_blank" rel="noreferrer">source</a></p>`,7))]),s("details",d,[s("summary",null,[i[3]||(i[3]=s("a",{id:"UncertaintyQuantification.SubSetInfinity",href:"#UncertaintyQuantification.SubSetInfinity"},[s("span",{class:"jlbinding"},"UncertaintyQuantification.SubSetInfinity")],-1)),i[4]||(i[4]=e()),n(t,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[5]||(i[5]=a(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SubSetInfinity</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(n</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, target</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Float64</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, levels</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, s</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Real</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Defines the properties of a Subset-∞ simulation where <code>n</code> is the number of initial samples, <code>target</code> is the target probability of failure at each level, <code>levels</code> is the maximum number of levels and <code>s</code> is the standard deviation for the proposal samples.</p><p><strong>Examples</strong></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> SubSetInfinity</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SubSetInfinity</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">100</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><strong>References</strong></p><p>[<a href="/UncertaintyQuantification.jl/previews/PR243/references#auRareEventSimulation2016">15</a>]</p><p>[<a href="/UncertaintyQuantification.jl/previews/PR243/references#patelliEfficientMonteCarlo2015">22</a>]</p><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/9573d7c5373e6153b9cbb7558b521a09977dfdb7/src/simulations/subset.jl#L43-L62" target="_blank" rel="noreferrer">source</a></p>`,8))]),s("details",g,[s("summary",null,[i[6]||(i[6]=s("a",{id:"UncertaintyQuantification.SubSetInfinityAdaptive",href:"#UncertaintyQuantification.SubSetInfinityAdaptive"},[s("span",{class:"jlbinding"},"UncertaintyQuantification.SubSetInfinityAdaptive")],-1)),i[7]||(i[7]=e()),n(t,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[8]||(i[8]=a(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SubSetInfinityAdaptive</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(n</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, target</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Float64</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, levels</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, Na</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, λ</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Real</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, s</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Real</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Implementation of: Papaioannou, Iason, et al. &quot;MCMC algorithms for subset simulation.&quot; Probabilistic Engineering Mechanics 41 (2015): 89-103</p><p>Defines the properties of a Subset-∞ adaptive where <code>n</code> are the number of samples per level, <code>target</code> is the target probability of failure at each level, <code>levels</code> is the maximum number of levels, <code>λ</code> (λ = 1 recommended) is the initial scaling parameter, and <code>Na</code> is the number simulations that will be run before <code>λ</code> is updated. Note that Na must be a multiple of n * target: <code>mod(ceil(n * target), Na) == 0)</code>. The initial variance of the proposal distribution is <code>s</code>.</p><p>Idea behind this algorithm is to adaptively select the correlation parameter of <code>s</code> at each intermediate level, by simulating a subset Na of the chains (which must be choosen without replacement at random) and modifying the acceptance rate towards the optimal αstar = 0.44</p><p><strong>Constructors</strong></p><ul><li><p><code>SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer)</code> (default: λ = s = 1)</p></li><li><p><code>SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real)</code> (λ = s)</p></li><li><p><code>SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real, s::Real)</code></p></li></ul><p><strong>Note</strong></p><p>The following constructors will run the same number of samples, but SubSetInfinityAdaptive will update <code>s</code> after each chain:</p><ul><li><p><code>SubSetInfinityAdaptive(400, 0.1, 10, 40)</code></p></li><li><p><code>SubSetInfinity(400, 0.1, 10, 0.5)</code></p></li></ul><p><strong>Examples</strong></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> SubSetInfinityAdaptive</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">200</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SubSetInfinityAdaptive</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">200</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p><strong>References</strong></p><p>[<a href="/UncertaintyQuantification.jl/previews/PR243/references#papaioannou2015mcmc">23</a>]</p><p>[<a href="/UncertaintyQuantification.jl/previews/PR243/references#chan2022adaptive">24</a>]</p><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/9573d7c5373e6153b9cbb7558b521a09977dfdb7/src/simulations/subset.jl#L92-L132" target="_blank" rel="noreferrer">source</a></p>`,15))])])}const m=l(r,[["render",E]]);export{f as __pageData,m as default};
