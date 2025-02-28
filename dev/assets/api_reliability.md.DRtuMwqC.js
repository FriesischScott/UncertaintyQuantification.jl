import{_ as l,C as o,c as r,o as p,ai as a,j as t,a as s,G as n}from"./chunks/framework.CmMe_2x7.js";const m=JSON.parse('{"title":"Reliability","description":"","frontmatter":{},"headers":[],"relativePath":"api/reliability.md","filePath":"api/reliability.md","lastUpdated":null}'),d={name:"api/reliability.md"},h={class:"jldocstring custom-block",open:""},c={class:"jldocstring custom-block",open:""},k={class:"jldocstring custom-block",open:""};function u(y,i,f,b,g,F){const e=o("Badge");return p(),r("div",null,[i[9]||(i[9]=a('<h1 id="reliability" tabindex="-1">Reliability <a class="header-anchor" href="#reliability" aria-label="Permalink to &quot;Reliability&quot;">​</a></h1><h2 id="index" tabindex="-1">Index <a class="header-anchor" href="#index" aria-label="Permalink to &quot;Index&quot;">​</a></h2><ul><li><a href="#UncertaintyQuantification.FORM"><code>UncertaintyQuantification.FORM</code></a></li><li><a href="#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{&lt;:UQModel}}, Function, Union{UQInput, Vector{&lt;:UQInput}}, AbstractMonteCarlo}"><code>UncertaintyQuantification.probability_of_failure</code></a></li><li><a href="#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{&lt;:UQModel}}, Function, Union{UQInput, Vector{&lt;:UQInput}}, FORM}"><code>UncertaintyQuantification.probability_of_failure</code></a></li></ul><h2 id="types" tabindex="-1">Types <a class="header-anchor" href="#types" aria-label="Permalink to &quot;Types&quot;">​</a></h2>',4)),t("details",h,[t("summary",null,[i[0]||(i[0]=t("a",{id:"UncertaintyQuantification.FORM",href:"#UncertaintyQuantification.FORM"},[t("span",{class:"jlbinding"},"UncertaintyQuantification.FORM")],-1)),i[1]||(i[1]=s()),n(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),i[2]||(i[2]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">FORM</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(n</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">10</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,tol</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Real</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1e-3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,fdm</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">FiniteDifferencesMethod</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">CentralFiniteDifferences</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span></code></pre></div><p>used to perform the first order reliability method using the HLRF algorithm with <code>n</code> iterations and tolerance <code>tol</code>. Gradients are estimated through <code>fdm</code>.</p><p><strong>References</strong></p><p>[<a href="/UncertaintyQuantification.jl/dev/references#rackwitzStructuralReliability1978">12</a>]</p><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f5ee6cce729f0d6a57979257379c942cdf42f86f/src/reliability/form.jl#L1-L9" target="_blank" rel="noreferrer">source</a></p>',5))]),i[10]||(i[10]=t("h2",{id:"methods",tabindex:"-1"},[s("Methods "),t("a",{class:"header-anchor",href:"#methods","aria-label":'Permalink to "Methods"'},"​")],-1)),t("details",c,[t("summary",null,[i[3]||(i[3]=t("a",{id:"UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, FORM}",href:"#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, FORM}"},[t("span",{class:"jlbinding"},"UncertaintyQuantification.probability_of_failure")],-1)),i[4]||(i[4]=s()),n(e,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),i[5]||(i[5]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">probability_of_failure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(models</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Union{Vector{&lt;:UQModel},UQModel}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,performance</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">),inputs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Union{Vector{&lt;:UQInput},UQInput}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,sim</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">FORM</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Perform a reliability analysis using the first order reliability method (FORM), see <a href="/UncertaintyQuantification.jl/dev/api/reliability#UncertaintyQuantification.FORM"><code>FORM</code></a>. Returns the estimated probability of failure <code>pf</code>, the reliability index <code>β</code> and the design point <code>dp</code>.</p><p><strong>Examples</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>pf, β, dp = probability_of_failure(model, performance, inputs, sim)</span></span></code></pre></div><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f5ee6cce729f0d6a57979257379c942cdf42f86f/src/reliability/form.jl#L24-L34" target="_blank" rel="noreferrer">source</a></p>',5))]),t("details",k,[t("summary",null,[i[6]||(i[6]=t("a",{id:"UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, AbstractMonteCarlo}",href:"#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, AbstractMonteCarlo}"},[t("span",{class:"jlbinding"},"UncertaintyQuantification.probability_of_failure")],-1)),i[7]||(i[7]=s()),n(e,{type:"info",class:"jlObjectType jlMethod",text:"Method"})]),i[8]||(i[8]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">probability_of_failure</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(models</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Union{Vector{&lt;:UQModel},UQModel}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,performance</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">),inputs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Union{Vector{&lt;:UQInput},UQInput}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,sim</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">AbstractMonteCarlo</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Perform a reliability analysis with a standard Monte Carlo simulation. Returns the estimated probability of failure <code>pf</code>, the standard deviation <code>σ</code> and the <code>DataFrame</code> containing the evaluated <code>samples</code>. The simulation <code>sim</code> can be any instance of <code>AbstractMonteCarlo</code>.</p><p><strong>Examples</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>pf, σ, samples = probability_of_failure(model, performance, inputs, sim)</span></span></code></pre></div><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f5ee6cce729f0d6a57979257379c942cdf42f86f/src/reliability/probabilityoffailure.jl#L1-L12" target="_blank" rel="noreferrer">source</a></p>',5))])])}const U=l(d,[["render",u]]);export{m as __pageData,U as default};
