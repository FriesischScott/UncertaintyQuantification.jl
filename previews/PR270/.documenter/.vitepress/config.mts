import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/UncertaintyQuantification.jl/previews/PR270/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Home', link: '/index' },
{ text: 'Manual', collapsed: false, items: [
{ text: 'Introduction', link: '/manual/introduction' },
{ text: 'Getting Started', link: '/manual/gettingstarted' },
{ text: 'Kernel Density Estimation', link: '/manual/kde' },
{ text: 'Gaussian Mixture Models', link: '/manual/gaussianmixture' },
{ text: 'Reliability Analysis', link: '/manual/reliability' },
{ text: 'Metamodelling', link: '/manual/metamodels' },
{ text: 'Simulations', link: '/manual/simulations' },
{ text: 'Bayesian Updating', link: '/manual/bayesianupdating' },
{ text: 'Parallelisation', link: '/manual/parallelisation' },
{ text: 'Stochastic Dynamics', link: '/manual/dynamics' },
{ text: 'High Performance Computing', link: '/manual/hpc' }]
 },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Gaussian Mixture Model', link: '/examples/inputs' },
{ text: 'External Models', link: '/examples/external' },
{ text: 'Metamodels', link: '/examples/metamodels' },
{ text: 'Bayesian Updating', link: '/examples/bayesianupdating' },
{ text: 'Stochastic Dynamics', link: '/examples/dynamics' },
{ text: 'High Performance Computing', link: '/examples/hpc' }]
 },
{ text: 'Benchmarks', collapsed: false, items: [
{ text: 'Subset Simulation', link: '/benchmarks/subset' }]
 },
{ text: 'API', collapsed: false, items: [
{ text: 'Inputs', link: '/api/inputs' },
{ text: 'Models', link: '/api/models' },
{ text: 'Reliability', link: '/api/reliability' },
{ text: 'ResponseSurface', link: '/api/responsesurface' },
{ text: 'PolyharmonicSpline', link: '/api/polyharmonicspline' },
{ text: 'Simulations', link: '/api/simulations' },
{ text: 'Bayesian Updating', link: '/api/bayesianupdating' },
{ text: 'Power Spectral Density Functions', link: '/api/psd' },
{ text: 'Stochastic Processes (Spectral Representation)', link: '/api/spectralrepresentation' },
{ text: 'SlurmInterface', link: '/api/slurm' },
{ text: 'Polynomial Chaos Expansions', link: '/api/polynomialchaos' }]
 },
{ text: 'References', link: '/references' }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker',
  }
]
// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/UncertaintyQuantification.jl/previews/PR270/', // TODO: replace this in makedocs!
  title: 'UncertaintyQuantification.jl',
  description: 'Documentation for UncertaintyQuantification.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [
    
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  vite: {
    build: {
      assetsInlineLimit: 0, // so we can tell whether we have created inlined images or not, we don't let vite inline them
    }
  },

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  themeConfig: {
    outline: 'deep',
    // https://vitepress.dev/reference/default-theme-config
    
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Manual', collapsed: false, items: [
{ text: 'Introduction', link: '/manual/introduction' },
{ text: 'Getting Started', link: '/manual/gettingstarted' },
{ text: 'Kernel Density Estimation', link: '/manual/kde' },
{ text: 'Gaussian Mixture Models', link: '/manual/gaussianmixture' },
{ text: 'Reliability Analysis', link: '/manual/reliability' },
{ text: 'Metamodelling', link: '/manual/metamodels' },
{ text: 'Simulations', link: '/manual/simulations' },
{ text: 'Bayesian Updating', link: '/manual/bayesianupdating' },
{ text: 'Parallelisation', link: '/manual/parallelisation' },
{ text: 'Stochastic Dynamics', link: '/manual/dynamics' },
{ text: 'High Performance Computing', link: '/manual/hpc' }]
 },
{ text: 'Examples', collapsed: false, items: [
{ text: 'Gaussian Mixture Model', link: '/examples/inputs' },
{ text: 'External Models', link: '/examples/external' },
{ text: 'Metamodels', link: '/examples/metamodels' },
{ text: 'Bayesian Updating', link: '/examples/bayesianupdating' },
{ text: 'Stochastic Dynamics', link: '/examples/dynamics' },
{ text: 'High Performance Computing', link: '/examples/hpc' }]
 },
{ text: 'Benchmarks', collapsed: false, items: [
{ text: 'Subset Simulation', link: '/benchmarks/subset' }]
 },
{ text: 'API', collapsed: false, items: [
{ text: 'Inputs', link: '/api/inputs' },
{ text: 'Models', link: '/api/models' },
{ text: 'Reliability', link: '/api/reliability' },
{ text: 'ResponseSurface', link: '/api/responsesurface' },
{ text: 'PolyharmonicSpline', link: '/api/polyharmonicspline' },
{ text: 'Simulations', link: '/api/simulations' },
{ text: 'Bayesian Updating', link: '/api/bayesianupdating' },
{ text: 'Power Spectral Density Functions', link: '/api/psd' },
{ text: 'Stochastic Processes (Spectral Representation)', link: '/api/spectralrepresentation' },
{ text: 'SlurmInterface', link: '/api/slurm' },
{ text: 'Polynomial Chaos Expansions', link: '/api/polynomialchaos' }]
 },
{ text: 'References', link: '/references' }
]
,
    editLink: { pattern: "https://github.com/FriesischScott/UncertaintyQuantification.jl/edit/master/docs/src/:path" },
    socialLinks: [
      { icon: 'slack', link: 'https://julialang.org/slack/' }
    ],
    footer: {
      message: 'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> and <a href="https://luxdl.github.io/DocumenterVitepress.jl/stable/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
