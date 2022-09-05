<body><div class="wrapper">
	<header>
		<h1><a href="https://jerkwin.github.io/gmxtools/">gmxtools</a></h1>
		<p>scripts/programs/utilities/tools for GROMACS</p>
		<p class="view"><a href="#Introduction">Introduction</a></p>
		<p class="view"><a href="#Tools">Tools</a></p>
	</header>

<h2 id="Introduction">Introduction</h2>

<p>This is a collection of the tools I used for GROMACS. Hope it helps for you as well.</p>
<p>If you use one/some of the tools for your publication(s), I will appreciate if you cite the Zenodo DOI <a href="https://doi.org/10.5281/zenodo.6408973">10.5281/zenodo.6408973</a></p>

<h2 id="Tools">Tools</h2></p>

<table>
<th>名称<br>NAME</th><th>说明<br>INSTRUCTION</th> <th>代码<br>CODE</th> <th>示例<br>EXAMPLE</th>

<tr>
<td>calc</td>
<td><a href="https://jerkwin.github.io/2019/11/19/%E5%9C%A8%E7%BA%BF%E8%AE%A1%E7%AE%97%E5%85%BC%E5%8D%95%E4%BD%8D%E6%8D%A2%E7%AE%97%E5%99%A8/">在线计算兼单位换算器</a></td>
<td><a href="./calc/calc.html">calc</a><br>
<a href="./calc/gnuecharts.html">gnuECharts</a></td>
<td>N/A</td>
</tr>

<tr>
<td>wvmd</td>
<td><a href="https://jerkwin.github.io/2017/02/11/%E8%87%AA%E5%8A%A8%E8%B0%83%E6%95%B4VMD%E7%AA%97%E5%8F%A3%E7%9A%84%E4%BD%8D%E7%BD%AE%E5%92%8C%E5%A4%A7%E5%B0%8F/">自动调整VMD窗口的位置和大小</a></td>
<td><a href="./wvmd/wvmd.zip">wvmd.zip</a></td>
<td>N/A</td>
</tr>

<tr>
<td>graphene</td>
<td><a href="http://jerkwin.github.io/2014/05/09/%E7%9F%B3%E5%A2%A8%E7%83%AF-%E5%BB%BA%E6%A8%A1-%E5%87%A0%E4%BD%95%E6%80%A7%E8%B4%A8%E5%8F%8A%E5%8A%9B%E5%9C%BA%E6%A8%A1%E6%8B%9F/">石墨烯：建模, 几何性质及力场模拟</a><br>
<a href="https://jerkwin.github.io/2014/12/24/%E7%9F%B3%E5%A2%A8%E7%83%AF%E5%9C%A8%E7%BA%BF%E5%88%9B%E5%BB%BA%E5%B7%A5%E5%85%B7/">石墨烯在线创建工具</a><br>
<a href="https://jerkwin.github.io/GMX/GMXtut-8/">创建周期性体系的拓扑文件：以石墨烯为例</a><br>
<a href="https://jerkwin.github.io/2020/04/05/%E6%B0%A7%E5%8C%96%E7%9F%B3%E5%A2%A8%E7%83%AF%E7%9A%84%E7%BB%93%E6%9E%84%E4%B8%8E%E5%BB%BA%E6%A8%A1/">氧化石墨烯的结构与建模</a><br>
<a href="https://jerkwin.github.io/2020/06/03/%E7%9F%B3%E5%A2%A8%E7%83%AF%E5%B9%B3%E9%9D%A2%E7%9A%84%E5%8D%B7%E6%9B%B2%E4%B8%8E%E6%89%AD%E6%9B%B2/">石墨烯平面的卷曲与扭曲</a><br>
</td>
<td><a href="./model/graphene.html">graphene</a></td>
<td>N/A</td>
</tr>

<tr>
<td>pipistack</td>
<td><a href="https://jerkwin.github.io/2018/08/29/Pi-Pi%E5%A0%86%E7%A7%AF%E8%B7%9D%E7%A6%BB%E5%92%8C%E5%A0%86%E7%A7%AF%E8%A7%92%E5%BA%A6%E7%9A%84%E8%AE%A1%E7%AE%97/">Pi-Pi堆积距离和堆积角度的计算</a></td>
<td><a href="./pipistack/pipistack.tcl">pipistack.tcl</a><br>
<a href="./pipistack/pipistack_linalg.tcl">pipistack_linalg.tcl</a></td>
<td><a href="./pipistack/ph2.zip">ph2.zip</a></td>

</td>

<tr>
<td>fitmol</td>
<td><a href="./fitmol/README.md.html">README</a></td>
<td><a href="./fitmol/FitMol.f90">FitMol.f90</a></td>
<td><a href="./fitmol/ch4.xyz">ch4.xyz</a><br>
<a href="./fitmol/ch4_ch3oh.xyz">ch4_ch3oh.xyz</a>
</td>

<tr>
<td>xpm2all</td>
<td><a href="https://jerkwin.github.io/2018/05/09/xpm%E6%96%87%E4%BB%B6%E5%A4%84%E7%90%86%E8%84%9A%E6%9C%AC/">xpm文件处理脚本</a><br>
<a href="https://jerkwin.github.io/2020/02/29/%E5%88%86%E5%AD%90%E6%A8%A1%E6%8B%9F%E5%91%A8%E5%88%8A-%E7%AC%AC_8_%E6%9C%9F/">xpm2all: xpm文件转换工具</a><br>
<a href="http://jerkwin.github.io/2020/07/10/%E4%BD%BF%E7%94%A8xpm2all%E8%84%9A%E6%9C%AC%E8%AE%A1%E7%AE%97%E8%9B%8B%E7%99%BD%E4%BA%8C%E7%BA%A7%E7%BB%93%E6%9E%84%E6%BC%94%E5%8C%96%E5%8F%8A%E5%90%AB%E9%87%8F/">使用xpm2all脚本计算蛋白二级结构演化及含量</a><br>
<a href="https://jerkwin.github.io/2021/03/30/xpm2all%E6%9B%B4%E6%96%B0-%E4%BA%8C%E7%BA%A7%E7%BB%93%E6%9E%84%E7%BB%98%E5%88%B6,_%E9%A2%9C%E8%89%B2%E6%96%B9%E6%A1%88/">xpm2all更新：二级结构绘制, 颜色方案</a>

</td>
<td><a href="./xpm2all/xpm2all.bsh">xpm2all.bsh</a><br>
<a href="./xpm2all/xpm2rgb.bat">xpm2rgb.bat</a>
</td>
<td><a href="./xpm2all/ss.xpm">ss.xpm</a><br>
<a href="./xpm2all/gibbs.xpm">gibbs.xpm</a>
</td>

<tr>
<td>gmx_mmpbsa</td>
<td>
<a href="https://jerkwin.github.io/2019/07/31/gmx_mmpbsa%E4%BD%BF%E7%94%A8%E8%AF%B4%E6%98%8E/">gmx_mmpbsa使用说明</a><br>
<a href="https://jerkwin.github.io/2021/03/16/gmx_mmpbsa%E8%84%9A%E6%9C%AC%E6%9B%B4%E6%96%B0-%E5%B1%8F%E8%94%BD%E6%95%88%E5%BA%94%E4%B8%8E%E7%86%B5%E8%B4%A1%E7%8C%AE/">gmx_mmpbsa脚本更新：屏蔽效应与熵贡献</a><br>
<a href="https://jerkwin.github.io/2021/11/26/gmx_mmpbsa%E8%84%9A%E6%9C%AC%E6%9B%B4%E6%96%B0-%E6%B8%85%E7%90%86%E6%95%B4%E7%90%86%E8%BE%93%E5%87%BA/">gmx_mmpbsa脚本更新：清理整理输出</a><br>
<a href="http://jerkwin.github.io/2022/02/09/gmx_mmpbsa%E6%9B%B4%E6%96%B0-%E4%BF%AE%E6%AD%A3%E5%8E%9F%E5%AD%90%E5%8D%8A%E5%BE%84bug,_%E6%94%B9%E7%94%A8AMBRR_PB4,_%E4%B8%99%E6%B0%A8%E9%85%B8%E6%89%AB%E6%8F%8FCAS/">gmx_mmpbsa更新：修正原子半径bug, 改用AMBRR PB4, 丙氨酸扫描CAS</a>

</td>
<td><a href="./gmx_mmpbsa/gmx_mmpbsa.bsh">gmx_mmpbsa.bsh</a></td>
<td><a href="./gmx_mmpbsa/1ebz.zip">1ebz.zip</a></td>
</tr>

<tr>
<td>gmx_ir</td>
<td><a href="https://jerkwin.github.io/2017/08/20/%E4%BD%BF%E7%94%A8GROMACS%E8%AE%A1%E7%AE%97%E7%BA%A2%E5%A4%96%E5%85%89%E8%B0%B1/">使用GROMACS计算红外光谱</a></td>
<td><a href="./gmx_ir/gmx_ir.bsh">gmx_ir.bsh</a></td>
<td>N/A</td>
</tr>

<tr>
<td>gmx_hbdat</td>
<td>
<a href="https://jerkwin.github.io/2021/06/19/GROMACS氢键分析工具hbond的使用及扩展/">GROMACS氢键分析工具hbond的使用及扩展</a><br>
<a href="https://jerkwin.github.io/2022/05/23/使用gnuplot绘制氢键数据的平行轴图/">使用gnuplot绘制氢键数据的平行轴图</a><br>
</td>
<td><a href="./gmx_hbdat/gmx_hbdat.bsh">gmx_hbdat.bsh</a></td>
<td><a href="./gmx_hbdat/hbdat_1crn.zip">hbdat_1crn.zip</a></td>
</tr>

<tr>
<td>remd_tgenerator</td>
<td><a href="https://jerkwin.github.io/2021/09/30/副本交换动力学T-REMD模拟的温度分布计算器/">副本交换动力学T-REMD模拟的温度分布计算器</a></td>
<td><a href="./remd_tgenerator/remd_tgenerator.html">remd_tgenerator</a></td>
<td>N/A</td>
</tr>

<tr>
<td>xff</td>
<td><a href="http://jerkwin.github.io/2022/04/15/%E5%8A%9B%E5%9C%BA%E6%8B%9F%E5%90%88%E5%B7%A5%E5%85%B7xff%E5%BC%80%E5%8F%91%E6%9D%82%E8%AE%B0/">力场拟合工具xff开发杂记</a></td>
<td><a href="./xff/xff.html">xff</a></td>
<td><a href="./xff/c6h6.fchk">c6h6.fchk</a><br>
<a href="./xff/sf6.FChk">sf6.FChk</a>
</td>
</tr>

<tr>
<td>cagen</td>
<td><a href="https://jerkwin.github.io/2022/09/02/%E6%B0%B4%E6%BA%B6%E6%B6%B2%E4%BD%93%E7%B3%BB%E6%88%90%E7%AC%BC%E5%88%86%E6%9E%90%E7%A8%8B%E5%BA%8Fcagen/">水溶液体系成笼分析程序cagen</a></td>
<td>N/A</td>
<td><a href="./cagen/cagen.zip">cagen.zip</a></td>
</tr>

</table>

<h2 id="Questions">Questions?</h2></p>

If you have any quesiotns about these tools, please let me know.

<ul>
<li><a href="https://groups.google.com/forum/#!forum/gmxtools">Google Forum</a></li>
<li><a href="https://github.com/Jerkwin/gmxtools/issues">GitHub issues</a></li>
</ul>

</div></body>