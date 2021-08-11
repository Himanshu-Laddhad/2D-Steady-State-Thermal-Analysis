<h1>2D Plate Thermal analysis </h1>
<h2>Introduction</h2>
<p>
  This project calculates the temperature distribution of cell centers of a 2D plate. This is done by solving the steady state equation for the given 2D domain by finite volume method. Uniform and Non-Uniform meshes in a increasing/decreasing/constant cell length fashion can be created by putting suitable value of successive ratio.
</p> 
<break>
<h2>Requirements</h2>
<p>
  The basic information that has to be provided to the program include:
  <ul>
    <li>Length</li>
    <li>Width</li>
    <li>Successive ratio in x direction</li>
    <li>Successive ratio in y direction</li>
    <li>Cells in x direction</li>
    <li>Cells in y direction</li>
    <li>Boundary coditions</li>
    <li>Number of iterations</li>
  </ul>
 For graph plotting, the information required is as follows:
  <ul>
    <li>Contour plot is plotted automatically.</li>
    <li>3d Surface plot is plotted automatically.</li>
    <li>Temperature distribution for particular value of x - cell number (in x direction) is required</li>
    <li>Temperature distribution for particular value of y - cell number (in y direction) is required</li>
    <li>Cells in x direction</li>
  </ul>
</p>
<h2>Example</h2>
<h3>Meshes:</h3>
<p>Different kind of meshes that can be generated are as follows: </p>
<img src = "https://user-images.githubusercontent.com/63182419/128965926-6bd5b95a-e644-4b33-89ab-32e97afe12f4.png"></img>
<p>Uniform Square Mesh</p>
<img src = "https://user-images.githubusercontent.com/63182419/128965928-71a0f28c-cdd0-4dc9-a208-06bfa00fbc2e.png"></img>
<p>Uniform Rectangular Mesh</p>
<img src = "https://user-images.githubusercontent.com/63182419/128965930-84635c67-54de-4f9d-bf29-3cf066a3917d.png"></img>
<p>Non-Uniform Rectangular Mesh</p>
<break>
<h3>Graphs:</h3>
  <p>Data Input
    <ul>
    <li>Length: 40</li>
    <li>Width: 40</li>
    <li>Successive ratio in x direction: 1</li>
    <li>Successive ratio in y direction: 1</li>
    <li>Cells in x direction: 40</li>
    <li>Cells in y direction: 40</li>
    <li>Boundary coditions:
      <ul>
        <li>Left boundary: 100</li>
        <li>Right boundary: 0</li>
        <li>Top boundary: 0</li>
        <li>Bottom boundary: 0</li>
    </li>
    <li>Number of iterations: 400</li>
  </ul>
 </p>
<img src = "https://user-images.githubusercontent.com/63182419/128967032-d75af0f9-4c16-4be0-a09e-ca5b9243667d.png"></img>
<p>3D Temperature Surface plot</p>
<img src = "https://user-images.githubusercontent.com/63182419/128967035-f2eab516-910c-40e1-bc22-1524555688e2.png"></img>
<p>Temperature Contour</p>
<img src = "https://user-images.githubusercontent.com/63182419/128967039-e93f1d7c-2e9a-4822-9bdf-ffe4f79873b3.png"></img>
<p>Variation of Temp in x direc for 10th cell in 40cells</p>
<img src = "https://user-images.githubusercontent.com/63182419/128967027-f22c233c-6616-40d0-a1ae-c96e18061b39.png"></img>
<p>Variation of Temp in y direc for 10th cell in 40cells</p>

<footer> Note: 3D surface plot and Temp contour are not functional for rectangular meshes but they can be solved and the results are stored in TempM.</footer>



