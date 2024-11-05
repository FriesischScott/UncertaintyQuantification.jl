import{_ as e,c as n,a5 as a,j as i,a as l,G as r,B as h,o}from"./chunks/framework.DDSZQJhB.js";const f=JSON.parse('{"title":"SlurmInterface","description":"","frontmatter":{},"headers":[],"relativePath":"api/slurm.md","filePath":"api/slurm.md","lastUpdated":null}'),p={name:"api/slurm.md"},d={class:"jldocstring custom-block",open:""};function k(u,s,c,g,y,m){const t=h("Badge");return o(),n("div",null,[s[3]||(s[3]=a('<h1 id="slurminterface" tabindex="-1">SlurmInterface <a class="header-anchor" href="#slurminterface" aria-label="Permalink to &quot;SlurmInterface&quot;">​</a></h1><h2 id="index" tabindex="-1">Index <a class="header-anchor" href="#index" aria-label="Permalink to &quot;Index&quot;">​</a></h2><ul><li><a href="#UncertaintyQuantification.SlurmInterface"><code>UncertaintyQuantification.SlurmInterface</code></a></li></ul><h2 id="type" tabindex="-1">Type <a class="header-anchor" href="#type" aria-label="Permalink to &quot;Type&quot;">​</a></h2>',4)),i("details",d,[i("summary",null,[s[0]||(s[0]=i("a",{id:"UncertaintyQuantification.SlurmInterface",href:"#UncertaintyQuantification.SlurmInterface"},[i("span",{class:"jlbinding"},"UncertaintyQuantification.SlurmInterface")],-1)),s[1]||(s[1]=l()),r(t,{type:"info",class:"jlObjectType jlType",text:"Type"})]),s[2]||(s[2]=a(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SlurmInterface</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(options</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Dict{String,String}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, throttle</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, btachsize</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Integer</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, extras</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Vector{String}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>When <code>SlurmInterface</code> is passed to an <code>ExternalModel</code>, model evaluations are executed using slurm job arrays. This allows for heavier simulations or workflows to be sampled, without relying on Julia&#39;s native parallelism. <code>SlurmInterface</code> automatically generates a slurm job array script, and Julia waits for this job to finish before extracting results.</p><p>When using <code>SlurmInterface</code>, you no longer need to load workers into Julia with <code>addprocs(N)</code>, and the requested nodes / tasks those required by individual model evaluations. Use <code>extras</code> to specify anything that must be preloaded for your models to be executed (for example loading modules).</p><p>The <code>throttle</code> specifies the number of simulations in the job array which are run concurrently. I.e., if you perform <code>MonteCarlo</code> simulation with <code>N=1000</code> samples, with <code>throttle=200</code>, it will run 1000 simulations in total, but only 200 at the same time. Your HPC scheduler (and admin) may be unhappy if you request too many concurrent jobs. If left empty, you scheduler&#39;s default throttle will be used. In the case that your HPC machine limits the size of submitted job arrays, you can split the submissions into smaller &quot;batches&quot;. Specify &quot;batchsize&quot; to the maximum size of a job array. This does not change the total number of runs.</p><p><strong>parameters</strong></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>options   : A dictionary of SBATCH options to add to the slurm script</span></span>
<span class="line"><span>throttle  : the number of jobs to be run at the same time</span></span>
<span class="line"><span>batchsize : maximum size of the slurm array, use when HPC limits the number of jobs in arrays</span></span>
<span class="line"><span>extras    : instructions to be executed before the model is run, e.g. activating a python environment or loading modules</span></span></code></pre></div><p><strong>Examples</strong></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> slurm </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> SlurmInterface</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Dict</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;account&quot;</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =&gt;</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;HPC_account_1&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;partition&quot;</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =&gt;</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;CPU_partition&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), extras </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;load python3&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SlurmInterface</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Dict</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;account&quot;</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =&gt;</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;HPC_account_1&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;partition&quot;</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =&gt;</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;CPU_partition&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, [</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;load python3&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span></code></pre></div><p><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/2a71f67bff400fb58a3eb13ec6aaa58903eec706/src/hpc/slurm.jl#L1-L24" target="_blank" rel="noreferrer">source</a></p>`,9))])])}const F=e(p,[["render",k]]);export{f as __pageData,F as default};
