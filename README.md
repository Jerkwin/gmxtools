<body><div class="wrapper">
	<header>
		<h1><a href="https://jerkwin.github.io/filamentcn/">Filament中文网</a></h1>
		<p>收集/整理/翻译 Filament 资料/文档/教程</p>
		<p class="view"><a href="#缘起">缘起</a><small> 为什么</small>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
		<a href="#概述">概述</a><small> 一句话介绍</small></p>
		<p class="view"><a href="#文档">文档</a><small> 构建, 材质</small>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
		<a href="#说明">说明</a><small> API, PPT</small></p>
		<p class="view"><a href="#教程">教程</a><small> 三角, 球, 模型</small>&nbsp;&nbsp;&nbsp;&nbsp;
		<a href="#示例">示例</a><small> 贴图, 材质, 动态</small></p>
		<p><img src="./jscss/二维码.png" alt=""></p>
	</header>

<h2 id="缘起">缘起</h2>

<p>一直以来, 我都对计算机图形学很感兴趣, 但无奈与平时工作相距甚远, 外兼视力不佳, 所以只能关注一下, 没有深入研究.</p>
<p>目前我已经翻译完成了Filament最重要的两篇材质文档, 外加三篇教程</p>
<p>其中的两篇材质文档是我目前看过的 <strong>最好的</strong> PBR资料, 对PBR渲染从原理到实现细节都有详细说明. 根据这些说明, 再加上一些图形学的知识, 完全可以实现自己的PBR.</p>
<p>我毕竟只是业余学习过计算机图形学, 很多名词的翻译不一定符合业内习惯. 如果你发现有不符合习惯的地方, 请指出修正. 如果你发现有错误和不合理的地方, 更要指出. 谢谢.</p>
<p>此外, 我还建立了一个QQ群, 用以方便大家交流.</p>

<h2 id="概述">Filament概述</h2></p>

<p><a href="https://github.com/google/filament">Filament</a>是一个用C++编写的基于物理的实时渲染器. 它优先考虑移动平台, 但也可用于多个平台.</p>
<p><iframe src="https://google.github.io/filament/webgl/demo_suzanne.html"></iframe></p>
<p>我们尽量保持Filament体积小, 加载快, 并专注渲染的特性. 例如, Filament不会在运行时编译材质. 相反, 我们提供了一个命令行工具<a href="https://github.com/google/filament/tree/master/tools/matc"><code>matc</code></a>, 用于离线编译材质.</p>

<h2 id="文档">Filament文档</h2>

<ul class="incremental">
<li><a href="">下载构建</a></li>
<li><a href="Filament.md.html">材质设计</a> <a href="https://google.github.io/filament/Filament.md.html">(原始文档)</a></li>
<li><a href="filamentcn//Materials.md.html">材质概览</a> <a href="https://google.github.io/filament/Materials.md.html">(原始文档)</a></li>
<li><a href="Material_Properties.pdf">材质参考页</a>
</li>
</ul>

<h2 id="说明">Filament说明</h2>

<ul class="incremental">
<li><a href="">JavaScript API参考</a></li>
<li><a href="">WebGL Meetup讲演 (2018)</a></li>
</li>
</ul>

<h2 id="教程">Filament教程</h2>
<ul class="incremental">
<li><a href="">三角形</a></li>
<li><a href="">红色球</a></li>
<li><a href="">Blender猴头</a></li>
</ul>

<h2 id="示例">Filament示例</h2>
<ul class="incremental">
<li><a href="">贴图材质球</a></li>
<li><a href="">Blender猴头</a></li>
<li><a href="">基本纽结</a></li>
</ul>

</div></body>