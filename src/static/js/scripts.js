/*!
* Start Bootstrap - Bare v5.0.9 (https://startbootstrap.com/template/bare)
* Copyright 2013-2023 Start Bootstrap
* Licensed under MIT (https://github.com/StartBootstrap/startbootstrap-bare/blob/master/LICENSE)
*/
// This file is intentionally blank
// Use this file to add JavaScript to your project
document.addEventListener('DOMContentLoaded', function() {
  function toggleBox(id) {
    var content = document.getElementById(id);
    var computedStyle = window.getComputedStyle(content);

    if (computedStyle.display === 'none') {
    content.style.display = 'block';
      } else {
        content.style.display = 'none';
      }
  }

  var headers = document.getElementsByClassName('header');
      for (var i = 0; i < headers.length; i++) {
        headers[i].onclick = function() {
          toggleBox(this.getAttribute('data-target'));
        };
  }

  function validateDecimal(input) {
      // Get the value entered by the user
      var value = parseFloat(input.value);

      // Check if the value has more than two decimal places
      if (value && value % 0.01 !== 0) {
        // Round the value to two decimal places
        var roundedValue = Math.round(value * 100) / 100;

        // Update the input value with the rounded value
        input.value = roundedValue;
      }
  }

  // Function to validate integer input fields
  function validateInteger(input) {
      var value = parseInt(input.value);
      if (Number.isInteger(value)) {
        input.value = value;
      } else {
        var roundedValue = Math.round(value);
        input.value = roundedValue;
      }
    }

  // Attach the validateDecimal function to the number input fields with step=0.01
  var decimalInputs = document.querySelectorAll('input[type="number"][step="0.01"]');
  for (var j = 0; j < decimalInputs.length; j++) {
    decimalInputs[j].addEventListener('input', function() {
      validateDecimal(this);
    });
  }

  // Attach the validateInteger function to the number input fields with step=1
  var integerInputs = document.querySelectorAll('input[type="number"][step="1"]');
    for (var k = 0; k < integerInputs.length; k++) {
      integerInputs[k].addEventListener('input', function() {
        validateInteger(this);
      });
    }

  // Get the input elements
  var inputs = document.querySelectorAll('input[type="number"]');

  // Add event listeners to the input elements
  inputs.forEach(function(input) {
    var step = parseFloat(input.getAttribute('step'));
    if (step % 1 === 0) {
      // For integer steps, bind the validateInteger function
      input.addEventListener('input', validateInteger);
    } else {
      // For decimal steps, bind the validateDecimal function
      input.addEventListener('input', validateDecimal);
    }
  });


  function checkNPEntryCount() {
  var table = document.getElementById('data-table-NP');
  var numRows = table.getElementsByTagName('tbody')[0].getElementsByTagName('tr').length;
  var submitNPButton = document.getElementsByName('submit_NP')[0]; // Update this line

  if (numRows > 2) {
    submitNPButton.disabled = true;
  } else {
    submitNPButton.disabled = false;
  }
  }

  function checkBREntryCount() {
  var table = document.getElementById('data-table-BR');
  var numRows = table.getElementsByTagName('tbody')[0].getElementsByTagName('tr').length;
  var submitBRButton = document.getElementsByName('submit_BR')[0]; // Update this line
  if (numRows > 4) {
    submitBRButton.disabled = true;
  } else {
    submitBRButton.disabled = false;
  }
  }


function cubesIntersect(p1, p2, p1_check, p2_check) {
    // Check for overlap in the X dimension
    if (p1[0] > p2_check[0] || p2[0] < p1_check[0]) {
        return true; // No overlap in the X dimension, so cubes do not intersect
    }

    // Check for overlap in the Y dimension
    if (p1[1] > p2_check[1] || p2[1] < p1_check[1]) {
        return true; // No overlap in the Y dimension, so cubes do not intersect
    }

    // Check for overlap in the Z dimension
    if (p1[2] > p2_check[2] || p2[2] < p1_check[2]) {
        return true; // No overlap in the Z dimension, so cubes do not intersect
    }

    // If there is no overlap in any dimension, the cubes do not intersect
    return false;
}

function canAddNewNP(npTable, p1, p2) {
    // Perform your calculations here
    var numRows = npTable.getElementsByTagName('tbody')[0].rows.length;
    if (numRows >= 1) {
        // Loop over all rows to verify information:
        for (var i = 0; i < numRows; i++) {
            var row = npTable.getElementsByTagName('tbody')[0].rows[i];
            var cell1Value = row.cells[0].textContent;
            var cell2Value = row.cells[1].textContent;
            var p1_check = cell1Value.split(',').map(parseFloat);
            var p2_check = cell2Value.split(',').map(parseFloat);

            // Check if the cube defined by P1 and P2 intersects the cube defiend by p1_check and p2_check:
            var doCubesIntersect = cubesIntersect(p1, p2, p1_check, p2_check);
            if (doCubesIntersect === false) {
                // The cubes intersect, so you can return or perform further actions.
                return false;
            }
        }
        // If we never return false then they do not intersect and we return true:
        return true
    } else {
        return true  // This means that the table is empty and we can always add an NP
    }
}


// Function to update the 3D plot
function update3DPlot() {
    // First we plot the Unit Cell bounds based on the values of maxX, maxY, and maxZ:
    const maxX = document.getElementById('maxX').value;
    const maxY = document.getElementById('maxY').value;
    const maxZ = document.getElementById('maxZ').value;
    const corners = [
          [0, 0, 0], [maxX, 0, 0],
          [maxX, 0, maxZ], [0, 0, maxZ],
          [0, maxY, 0], [maxX, maxY, 0],
          [maxX, maxY, maxZ], [0, maxY, maxZ],
        ];
    const edges = [
          [0, 1], [1, 2], [2, 3], [3, 0],
          [4, 5], [5, 6], [6, 7], [7, 4],
          [0, 4], [1, 5], [2, 6], [3, 7],
        ];
    const edgeTraces = edges.map(edge => {
      const x = [corners[edge[0]][0], corners[edge[1]][0]];
      const y = [corners[edge[0]][1], corners[edge[1]][1]];
      const z = [corners[edge[0]][2], corners[edge[1]][2]];

      return {
        type: 'scatter3d',
        mode: 'lines',
        x: x,
        y: y,
        z: z,
        line: {
          color: 'black', // Set line color to black
          width: 6,       // You can adjust line width as needed
          dash: 'dash',   // Set the line style to dashed
        },
        showlegend: false,  // Hide the legend entry for this trace
      };
    });

    const data = [
      ...edgeTraces,
      // Add other traces or data if needed
    ];

    // Next we loop over the NP Table to add any NPs the user has defined:
    var npTable = document.getElementById('data-table-NP');
    var numRows = npTable.getElementsByTagName('tbody')[0].rows.length;

    if (numRows >= 1){
        for (var i = 0; i < numRows; i++) {
                var row = npTable.getElementsByTagName('tbody')[0].rows[i];
                var cell1Value = row.cells[0].textContent.replace(/[^\d.,-]/g, '');
                var cell2Value = row.cells[1].textContent.replace(/[^\d.,-]/g, '');
                var p1 = cell1Value.split(',').map(parseFloat);
                var p2 = cell2Value.split(',').map(parseFloat);
                // Code here to plot a cube from p1_check to p2_check
                // Create the vertices of the rectangular prism
              const vertices = [
                p1,
                [p2[0], p1[1], p1[2]],
                [p2[0], p2[1], p1[2]],
                [p1[0], p2[1], p1[2]],
                [p1[0], p1[1], p2[2]],
                [p2[0], p1[1], p2[2]],
                [p2[0], p2[1], p2[2]],
                [p1[0], p2[1], p2[2]],
              ];

              // Define the vertices and faces for the Mesh3D trace
              const x = vertices.map(vertex => vertex[0]);
              const y = vertices.map(vertex => vertex[1]);
              const z = vertices.map(vertex => vertex[2]);

              const faces = [
                [0, 1, 2], // Front face (triangulated)
                [0, 2, 3], // Front face (triangulated)
                [4, 5, 6], // Back face (triangulated)
                [4, 6, 7], // Back face (triangulated)
                [0, 1, 5], // Bottom face (triangulated)
                [0, 5, 4], // Bottom face (triangulated)
                [2, 3, 7], // Top face (triangulated)
                [2, 7, 6], // Top face (triangulated)
                [0, 3, 7], // Left face (triangulated)
                [0, 7, 4], // Left face (triangulated)
                [1, 2, 6], // Right face (triangulated)
                [1, 6, 5], // Right face (triangulated)
              ];
              // Define the plotly trace:
              const trace = {
                type: 'mesh3d',
                x: x,
                y: y,
                z: z,
                i: faces.map(face => face[0]),
                j: faces.map(face => face[1]),
                k: faces.map(face => face[2]),
                facecolor: 'red',
                flatshading: true,
                lighting: {
                    ambient: 0.5, // Adjust ambient light
                    diffuse: 0.5, // Adjust diffuse light
                    specular: 0.5, // Adjust specular light
                  },
              };
              // Add the trace to the data array for the plot
              data.push(trace);
    }}


    var brTable = document.getElementById('data-table-BR');
    var rows = brTable.getElementsByTagName('tbody')[0].rows;
    var numRows2 = rows.length;
        for (var i = 0; i < numRows2; i++) {
                console.log(numRows2)
                console.log(i)
                var row = rows[i]; // Access the row directly from the rows collection
                var cell1Value = row.cells[0].textContent.replace(/[^\d.,-]/g, '');
                var cell2Value = row.cells[1].textContent.replace(/[^\d.,-]/g, '');
                var cell3Value = row.cells[2].textContent.replace(/[^\d.,-]/g, '');
                var cell4Value = row.cells[3].textContent.replace(/[^\d.,-]/g, '');
                if (cell2Value === '' && cell3Value === '' && cell4Value === '') {
                    // This is a binding POINT (i.e. scatter point) to add
                    var p1 = cell1Value.split(',').map(parseFloat);
                    var trace = {
                      x: [p1[0]], // x-coordinate
                      y: [p1[1]], // y-coordinate
                      z: [p1[2]], // z-coordinate
                      mode: 'markers', // Set the trace mode to markers
                      type: 'scatter3d',
                      marker: {
                        color: 'green', // Marker color
                        symbol: 'square', // Marker symbol (square)
                        size: 10 // Marker size
                      },
                      showlegend: false,
                    };

                }
                else if (cell3Value === '' && cell4Value === '') {
                    // This is a binding EDGE (i.e. line) to add
                    var p1 = cell1Value.split(',').map(parseFloat);
                    var p2 = cell2Value.split(',').map(parseFloat);
                    var trace = {
                      x: [p1[0], p2[0]], // x-coordinate
                      y: [p1[1], p2[1]], // y-coordinate
                      z: [p1[2], p2[2]], // z-coordinate
                      mode: 'markers+lines', // Set the trace mode to markers
                      type: 'scatter3d',
                      marker: {
                        color: 'green', // Marker color
                        symbol: 'square', // Marker symbol (square)
                        size: 10 // Marker size
                      },
                      line: {
                          color: 'black', // Set line color to black
                          width: 6,       // You can adjust line width as needed
                        },
                      showlegend: false,
                    };
                }
                else if (cell4Value === '') {
                    // This is a binding FACE (i.e. mesh3d) to add
                    var p1 = cell1Value.split(',').map(parseFloat);
                    var p2 = cell2Value.split(',').map(parseFloat);
                    var p3 = cell3Value.split(',').map(parseFloat);
                    var x = [p1[0], p2[0], p3[0]];
                    var y = [p1[1], p2[1], p3[1]];
                    var z = [p1[2], p2[2], p3[2]];

                    // Define the vertices for the mesh
                    var vertices = {
                      x: x,
                      y: y,
                      z: z,
                    };

                    // Create the mesh3d trace
                    var trace = {
                      type: 'mesh3d',
                      x: vertices.x,
                      y: vertices.y,
                      z: vertices.z,
                      flatshading: true, // Enable flat shading
                      color: 'green', // Set the color to green
                      showlegend: false,
                    };
                }
                else {
                    // This is a binding FACE (i.e. mesh3d) to add
                    var p1 = cell1Value.split(',').map(parseFloat);
                    var p2 = cell2Value.split(',').map(parseFloat);
                    var p3 = cell3Value.split(',').map(parseFloat);
                    var p4 = cell4Value.split(',').map(parseFloat);

                    var x = [p1[0], p2[0], p3[0], p4[0]]; // X-coordinates
                    var y = [p1[1], p2[1], p3[1], p4[1]]; // Y-coordinates
                    var z = [p1[2], p2[2], p3[2], p4[2]]; // Z-coordinates
                    var Ti = [0, 0, 1, 1];
                    var j = [1, 2, 2, 3];
                    var k = [2, 3, 3, 0];

                    var trace = {
                      "type": "mesh3d",
                      "x": x,
                      "y": y,
                      "z": z,
                      "i": Ti,
                      "j": j,
                      "k": k,
                      "flatshading": true,
                      "color": "green",
                      "showlegend": false
                    }
                }
                data.push(trace); // After the proper trace is created we add to data for plotting

    }

    const layout = {
      scene: {
        aspectmode: 'data',
        xaxis: {
          title: 'X (nm)',
          titlefont: {
            size: 16, // Set title font size
          },
          tickfont: {
            size: 14, // Set tick label font size to 12
          },
        },
        yaxis: {
          title: 'Y (nm)',
          titlefont: {
            size: 16,
          },
          tickfont: {
            size: 14,
          },
        },
        zaxis: {
          title: 'Z (nm)',
          titlefont: {
            size: 16,
          },
          tickfont: {
            size: 14,
          },
        },
      },
      title: {
        text: 'Start Criteria Preview',
        font: {
          size: 18, // Set title font size to 18
        },
      },
      lighting: {
        ambient: 0.9, // Adjust ambient light
        diffuse: 0.2, // Adjust diffuse light
        specular: 0.1, // Adjust specular light
      },
      width: 800, // Set the width of the plot
      height: 600, // Set the height of the plot
    };
    const config = {
      displayModeBar: false, // Disable interactivity, including saving images
    };

    var graphContainer = document.getElementById('graph-container');
    Plotly.newPlot(graphContainer, data, layout, config);
}

update3DPlot(); // Call this so that the page immediately renders


// Listener to update plot if X, Y, or Z are changed:
const maxXInput = document.getElementById("maxX");
const maxYInput = document.getElementById("maxY");
const maxZInput = document.getElementById("maxZ");

// Add an event listener to the "maxX" input
maxXInput.addEventListener("input", function() {
  // Call the "update3DPlot()" function whenever the input changes
  update3DPlot();
});
maxYInput.addEventListener("input", function() {
  // Call the "update3DPlot()" function whenever the input changes
  update3DPlot();
});
maxZInput.addEventListener("input", function() {
  // Call the "update3DPlot()" function whenever the input changes
  update3DPlot();
});


// Updating the NP Table:
document.getElementById('submit_NP').addEventListener('click', function() {
        const pointInput = document.getElementById('C1').value.replace(/[^\d.,-]/g, '');
        const pointInput2 = document.getElementById('C2').value.replace(/[^\d.,-]/g, '');
        const pointArray = pointInput.split(',').map(val => parseFloat(val.trim()));
        const pointArray2 = pointInput2.split(',').map(val2 => parseFloat(val2.trim()));

         if (
        pointArray.length === 3 && pointArray.every(Number.isFinite) &&
        pointArray2.length === 3 && pointArray2.every(Number.isFinite)
        ) {
            const [X1, Y1, Z1] = pointArray;
            const [X2, Y2, Z2] = pointArray2;
            // First we need to verify that the NP is within the acceptable ranges for X Y and Z:
            const maxX = document.getElementById('maxX').value;
            const maxY = document.getElementById('maxY').value;
            const maxZ = document.getElementById('maxZ').value;

            if (X1 < 0 || X1 > maxX || X2 < 0 || X2 > maxX) {
            alert('X1 and X2 must be within the acceptable range (0 to ' + maxX + ').');
            } else if (Y1 < 0 || Y1 > maxY || Y2 < 0 || Y2 > maxY) {
                alert('Y1 and Y2 must be within the acceptable range (0 to ' + maxY + ').');
            } else if (Z1 < 0 || Z1 > maxZ || Z2 < 0 || Z2 > maxZ) {
                alert('Z1 and Z2 must be within the acceptable range (0 to ' + maxZ + ').');
            } else if (Z1 > Z2 || Y1 > Y2 || X1 > Z2) {
                alert('C2 must be larger than C1!');
            } else {
                // Next, we need to verify that these NP regions are not intersecting cubes, and if so, we do not allow
                // this NP to be added!
                var npTable = document.getElementById('data-table-NP');

                var canAddNewRow = canAddNewNP(npTable, pointArray, pointArray2);
                if (canAddNewRow ){
                    // First we DISABLE the entry fields for maxX, maxY, and maxZ:
                    document.getElementById('maxX').disabled = true;
                    document.getElementById('maxY').disabled = true;
                    document.getElementById('maxZ').disabled = true;
                    document.querySelector('.disabled-text').style.display = 'inline'

                    // Next we add both pointArray and pointArray2 to the data-table-NP
                    // Get a reference to the table body
                    var tbody = document.querySelector('#data-table-NP tbody');

                    var newRow = tbody.insertRow(-1);  // -1 adds row to end of table
                    // Insert cell elements for the new row
                    var cell1 = newRow.insertCell(0);
                    var cell2 = newRow.insertCell(1);
                    var cell3 = newRow.insertCell(2);

                    // Set the values of the cells with X, Y, Z values
                    cell1.innerHTML = '(' + pointArray.join(', ') + ')'; // e.g., "X1, Y1, Z1"
                    cell2.innerHTML = '(' + pointArray2.join(', ') + ')'; // e.g., "X2, Y2, Z2"

                    // Create a button in the third cell to remove the row if needed
                    var removeButton = document.createElement('button');
                    removeButton.type = 'button';
                    removeButton.className = 'remove-np-button';
                    removeButton.innerHTML = 'Remove';
                    cell3.appendChild(removeButton);

                    // Lastly, after adding a NP, we call the update function:
                    checkNPEntryCount()  // Count how many binders are included to prevent overload
                    update3DPlot()
                }
                else {
                alert('Can not add this NP due to overlapping with an existing NP in the table!')
                }
                }

        } else {
            alert('Invalid input. Please enter valid 3D points for C1 and C2 as "X, Y, Z".');
        }
    });


// Listener to delete rows from the table if pressed
var table = document.getElementById('data-table-NP');
table.addEventListener('click', function(event) {
    if (event.target.classList.contains('remove-np-button')) {
        // Find the row that contains the clicked button
        var row = event.target.closest('tr');

        // Check if a valid row was found
        if (row) {
            // Delete the found row
            var tbody = table.getElementsByTagName('tbody')[0];
            var check_rows = row.rowIndex - 1
            tbody.deleteRow((row.rowIndex - 1));  // Delete the -1 index since the header row is there
            // Now if there are no more rows, then we re-enable the disabled
            var brTable = document.getElementById('data-table-BR');
            var brRowCount = brTable.getElementsByTagName('tbody')[0].rows.length;
            var npTable = document.getElementById('data-table-NP');
            var npRowCount = table.getElementsByTagName('tbody')[0].rows.length;
            if (npRowCount === 0 && brRowCount === 0) {
                // If both have 0 rows, then we re-enable the maxX, maxY, and maxZ values:
                document.querySelector('.disabled-text').style.display = 'none'
                document.getElementById('maxX').disabled = false;
                document.getElementById('maxY').disabled = false;
                document.getElementById('maxZ').disabled = false;
            }
            // Prevent further propagation of the click event
            checkNPEntryCount()  // Verify not too many NPs are included...
            update3DPlot()
        }
    }
});


function verifyBRPoints(p1, p2, p3, p4) {
    // Perform your calculations here
    var masterFlag = true;
    const maxXInput = document.getElementById("maxX");
    const maxYInput = document.getElementById("maxY");
    const maxZInput = document.getElementById("maxZ");
    if (p1.length === 3 && p1.every(Number.isFinite)) {
        // If the point is of the valid format, we need to make sure it is between the bounds of the problem:
        if (p1[0] < 0 || p1[0] > maxXInput || p1[1] < 0 || p1[1] > maxYInput || p1[2] < 0 || p1[2] > maxZInput) {
        return false; // If these conditions are not met, we return False meaning invalid
        }
        }

    else{
    return false
    }

    // Next verify p2 thru p4 (which are OPTIONAL):
    if (p2.length === 3 && p2.every(Number.isFinite) || isNaN(p2)) {
        if (p2[0] < 0 || p2[0] > maxXInput || p2[1] < 0 || p2[1] > maxYInput || p2[2] < 0 || p2[2] > maxZInput) {
        return false;
        }
        }
    else{
    return false
    }

    if (p3.length === 3 && p3.every(Number.isFinite) || isNaN(p3)) {
        if (p3[0] < 0 || p3[0] > maxXInput || p3[1] < 0 || p3[1] > maxYInput || p3[2] < 0 || p3[2] > maxZInput) {
        return false;
        }
        }
    else{
    return false
    }
    if (p4.length === 3 && p4.every(Number.isFinite) || isNaN(p4)) {
        if (p4[0] < 0 || p4[0] > maxXInput || p4[1] < 0 || p4[1] > maxYInput || p4[2] < 0 || p4[2] > maxZInput) {
        return false;
        }
        // In this case, we also need to verify that the points entered ARE CCW when input or else bad input

        }
    else{
    return false
    }
    return masterFlag
}

function invalidBROrder(p1, p2, p3, p4) {
// Next we need to verify that if the user entered in p1 and p2 and p4 as that is invalid (need to input sequentially)
    if ((p1.length === 3 && p3.length === 3) && (p2.length === 1)) {
    return false}
    if ((p1.length === 3 && p4.length === 3) && (p2.length === 1 || p3.length === 1)) {
    return false}
    return true
    }



// Updating the BR Table AJ ADD THE BR INTERSECTING CHECK CONDITION HERE
document.getElementById('submit_BR').addEventListener('click', function() {
        const pointInput1 = document.getElementById('pointInput1').value.replace(/[^\d.,-]/g, '');
        const pointInput2 = document.getElementById('pointInput2').value.replace(/[^\d.,-]/g, '');
        const pointInput3 = document.getElementById('pointInput3').value.replace(/[^\d.,-]/g, '');
        const pointInput4 = document.getElementById('pointInput4').value.replace(/[^\d.,-]/g, '');

        const pointArray1 = pointInput1.split(',').map(val1 => parseFloat(val1.trim()));
        const pointArray2 = pointInput2.split(',').map(val2 => parseFloat(val2.trim()));
        const pointArray3 = pointInput3.split(',').map(val3 => parseFloat(val3.trim()));
        const pointArray4 = pointInput4.split(',').map(val4 => parseFloat(val4.trim()));

        // First check and verify that the first point is valid since it is absolutely required:
        const check_input = verifyBRPoints(pointArray1, pointArray2, pointArray3, pointArray4)
        const check_BR_order = invalidBROrder(pointArray1, pointArray2, pointArray3, pointArray4)
        if (check_input && check_BR_order) {
                var npTable = document.getElementById('data-table-BR');
                //var canAddNewRow = canAddNewNP(npTable, pointArray, pointArray2);
                var canAddNewRow = true;  // Adding more tracking to the Ply_Graphs script because i suck at javascript
                if (canAddNewRow){
                    // First we DISABLE the entry fields for maxX, maxY, and maxZ:
                    document.getElementById('maxX').disabled = true;
                    document.getElementById('maxY').disabled = true;
                    document.getElementById('maxZ').disabled = true;
                    document.querySelector('.disabled-text').style.display = 'inline'

                    // Next we add both pointArray and pointArray2 to the data-table-NP
                    // Get a reference to the table body
                    var tbody = document.querySelector('#data-table-BR tbody');

                    var newRow = tbody.insertRow(-1);  // -1 adds row to end of table
                    // Insert cell elements for the new row
                    var cell1 = newRow.insertCell(0);
                    var cell2 = newRow.insertCell(1);
                    var cell3 = newRow.insertCell(2);
                    var cell4 = newRow.insertCell(3);
                    var cell5 = newRow.insertCell(4);

                    // Set the values of the cells with X, Y, Z values
                    cell1.innerHTML = '(' + pointArray1.join(', ') + ')'; // e.g., "X1, Y1, Z1"
                    if (pointArray2.length !== 3) { // If the length isn't 3 then it must be NaN
                        cell2.innerHTML = ''
                    }
                    else {
                        cell2.innerHTML = '(' + pointArray2.join(', ') + ')';
                    }

                    if (pointArray3.length !== 3) {
                        cell3.innerHTML = ''
                    }
                    else {
                        cell3.innerHTML = '(' + pointArray3.join(', ') + ')';
                    }

                    if (pointArray4.length !== 3) {
                        cell4.innerHTML = ''
                    }
                    else {
                        cell4.innerHTML = '(' + pointArray4.join(', ') + ')';
                    }

                    // Create a button in the third cell to remove the row if needed
                    var removeButton = document.createElement('button');
                    removeButton.type = 'button';
                    removeButton.className = 'remove-br-button';
                    removeButton.innerHTML = 'Remove';
                    cell5.appendChild(removeButton);

                    // Lastly, after adding a NP, we call the update function:
                    checkBREntryCount()  // Count how many binders are included to prevent overload
                    update3DPlot()
                }
                else {
                alert('Can not add this binding region due to overlapping with an existing NP in the table!')
                }
        }
        else if (check_BR_order === false) {
        alert('Invalid binding region entry. You must enter points in sequentially (i.e., you can NOT assign Point 1 and Point 3 without assigning Point 2');
        }
        else{
        alert('Invalid input for one of the binding region points. Please enter valid 3D point for at least Binding Point 1 in the form (X, Y, Z). Verify that no coordinate is outside the bounding box (dotted line) region.');
        }
    });


// Listener to delete rows from the table if pressed
var brTable = document.getElementById('data-table-BR');
brTable.addEventListener('click', function(event) {
    if (event.target.classList.contains('remove-br-button')) {
        // Find the row that contains the clicked button
        var row = event.target.closest('tr');

        // Check if a valid row was found
        if (row) {
            // Delete the found row
            var brTable = document.getElementById('data-table-BR');
            var tbody = brTable.getElementsByTagName('tbody')[0];
            var check_rows = row.rowIndex - 1
            tbody.deleteRow((row.rowIndex - 1));  // Delete the -1 index since the header row is there
            // Now if there are no more rows, then we re-enable the disabled
            var brRowCount = brTable.getElementsByTagName('tbody')[0].rows.length;
            var npTable = document.getElementById('data-table-NP');
            var npRowCount = table.getElementsByTagName('tbody')[0].rows.length;
            if (npRowCount === 0 && brRowCount === 0) {
                // If both have 0 rows, then we re-enable the maxX, maxY, and maxZ values:
                document.querySelector('.disabled-text').style.display = 'none'
                document.getElementById('maxX').disabled = false;
                document.getElementById('maxY').disabled = false;
                document.getElementById('maxZ').disabled = false;
            }
            // Prevent further propagation of the click event
            checkBREntryCount()  // Count how many binders are included to prevent overload
            update3DPlot()
        }
    }
});


function convertNPTableToJSON(){
    // This function will read the active NP table and return the string that is required for the JSON file input
    // Format: [((0, 15, 15), (20, 35, 35))] --> [((C1), (C2)), ((C1), (C2)), ...]
    var activeString = '[';  // The NP will always open with a [\
    // Now we read in row by row of the np table:
    var npTable = document.getElementById('data-table-NP');
    if (npTable) {
      // Get the table body (tbody) element
      var tbody = npTable.getElementsByTagName('tbody')[0];

      // Check if the tbody exists
      if (tbody) {
        // Get all the rows in the tbody
        var rows = tbody.getElementsByTagName('tr');

        // Loop through each row
        for (var i = 0; i < rows.length; i++) {
          var row = rows[i];

          // Get the cells (td elements) in the row
          var cells = row.getElementsByTagName('td');

          // Check if there are at least two cells in the row
          if (cells.length >= 2) {
            // Extract the text content from the first and second cells (Corner 1 and Corner 2)
            var corner1 = cells[0].textContent;
            var corner2 = cells[1].textContent;

            // Append the values to the activeString with a separator, for example, a comma
            activeString += '(' + corner1 + ',' + corner2 + '),';
          }
        }
      }
    }
    activeString += ']'
    // We return the active string with a closed bracket:
    return activeString
}


function convertBRTableToJSON(){
    // This function will read the active BR table and return the string that is required for the JSON file input
    // Format: [((25, 0, 25),None), ((15, 50, 25),(35, 50, 25)), ((35, 35, 50),(35, 15, 50),(15, 15, 50), (15, 35, 50)), ((25, 25, 0),None)]
    var activeString = '['
    // Now read in row by row of br table:
    var brTable = document.getElementById('data-table-BR');
    if (brTable) {
      // Get the table body (tbody) element
      var tbody = brTable.getElementsByTagName('tbody')[0];

      // Check if the tbody exists
      if (tbody) {
        // Get all the rows in the tbody
        var rows = tbody.getElementsByTagName('tr');

        // Loop through each row
        for (var i = 0; i < rows.length; i++) {
          var row = rows[i];

          // Get the cells (td elements) in the row
          var cells = row.getElementsByTagName('td');


          // Check if there are at least two cells in the row
          if (cells.length >= 4) {
            // Extract the text content from the first and second cells (Corner 1 and Corner 2)
            var V1 = cells[0].textContent;
            var V2 = cells[1].textContent;
            var V3 = cells[2].textContent;
            var V4 = cells[3].textContent;
            // Now, adding these to the string is a bit weird since it depends on the values of V2 - V4:
            if (V2 === ''){
                activeString += '(' + V1 + ',None),';
            }
            else if (V3 === '') {
                activeString += '(' + V1 + ',' + V2 + '),';
            }
            else if (V4 === '') {
                activeString += '(' + V1 + ',' + V2 + ',' + V3 + '),';
            }
            else {
                activeString += '(' + V1 + ',' + V2 + ',' + V3 + ',' + V4 + '),';
            }
          }
        }
      }
    }
    activeString += ']' // Close bracket to end string
    // We return the active string with a closed bracket:
    return activeString
}


// Add an event listener to the button
document.getElementById("download_input_file").addEventListener("click", function() {
    // Initial check to see if we can download a file
    var brTable = document.getElementById('data-table-BR');
    var brRowCount = brTable.getElementsByTagName('tbody')[0].rows.length;
    if (brRowCount === 0) {
        alert('Can not download input file, you must specify binding regions!')
    }

    // After the file is successfully downloaded
    // Access various information from the HTML page
    var batchNumber = "MOSAInputGenerator";
    var max_x = document.getElementById('maxX').value;
    var max_y = document.getElementById('maxY').value;
    var max_z = document.getElementById('maxZ').value;
    var NP_Regions = convertNPTableToJSON();
    var routing_software = 'DAEDALUS';
    var binding_regions = convertBRTableToJSON();
    var random_seed = document.getElementById('random_seed').value;
    var min_edge_length = document.getElementById('minEdgeLength').value;
    var extend_rule_distance = document.getElementById('extend_rule_distance').value;
    var max_number_bps = document.getElementById('maxBPs').value;
    var max_edge_multiplier = document.getElementById('max_edge_mult').value; //max_edge_mult
    var min_face_angle = document.getElementById('minAngle').value; //minAngle
    var C_repulsion = document.getElementById('repulsion_const').value; //repulsion_const
    var allow_self_intersecting_geometry = document.getElementById('allow_self_intersecting_geometries').value; //allow_self_intersecting_geometries
    var max_number_edges_at_node = document.getElementById('max_edges_at_node').value; //max_edges_at_node
    var repulsion_distance = document.getElementById('repulsion_distance').value; //repulsion_distance
    var cooling_schedule = document.getElementById('cooling_schedule').value; //cooling_schedule
    var acceptance_function = document.getElementById('acceptance_function').value; //acceptance_function
    var numDecimals = document.getElementById('numDecimals').value; //numDecimals
    var max_number_of_epochs = document.getElementById('max_epochs').value; //max_epochs
    var NT1 = document.getElementById('NT1').value; //NT1
    var NT2 = document.getElementById('NT2').value; //NT2
    var Na = document.getElementById('Na').value; //Na
    var r_b = document.getElementById('r_b').value; //r_b
    var N_Bi = document.getElementById('N_Bi').value; //N_Bi
    var N_Bi_LowerLimit = document.getElementById('N_Bi_LowerLim').value; //N_Bi_LowerLim
    var minimal_candidate_set_size = document.getElementById('min_set_size').value; //min_set_size
    var r_i = document.getElementById('r_i').value; //r_i
    var max_time_of_simulation_minutes = document.getElementById('MaxSimTime').value; //MaxSimTime
    var cooling_rate_geometric = document.getElementById('cooling_rate').value;//cooling_rate
    var delta_T = document.getElementById('triki_delta').value; //triki_delta
    var T_min = document.getElementById('min_temperature').value; //min_temperature

    // Create JSON Object:
    var jsonObject = {
        batchNumber,
        max_x,
        max_y,
        max_z,
        NP_Regions,
        routing_software,
        binding_regions,
        random_seed,
        min_edge_length,
        extend_rule_distance,
        max_number_bps,
        max_edge_multiplier,
        min_face_angle,
        C_repulsion,
        allow_self_intersecting_geometry,
        max_number_edges_at_node,
        repulsion_distance,
        cooling_schedule,
        acceptance_function,
        numDecimals,
        max_number_of_epochs,
        NT1,
        NT2,
        Na,
        r_b,
        N_Bi,
        N_Bi_LowerLimit,
        minimal_candidate_set_size,
        r_i,
        max_time_of_simulation_minutes,
        cooling_rate_geometric,
        delta_T,
        T_min
    };

    // Convert the JSON object to a JSON string
    var jsonData = JSON.stringify(jsonObject, null, 2); // The second argument is for pretty printing

    // Create a Blob and a download link
    var blob = new Blob([jsonData], { type: "application/json" });
    var url = URL.createObjectURL(blob);

    var a = document.createElement('a');
    a.href = url;
    a.download = 'MOSA_Input.json'; // You can specify the filename here
    a.style.display = 'none'; // Hide the download link
    a.click();
    URL.revokeObjectURL(url);
});


// Add event listener to preload example inputs:

function clearBRnNPTables(){
    // This function will clear out all rows inside of the NP and BR tables to make preloading examples easier:
    var npTable = document.getElementById("data-table-NP");
    var brTable = document.getElementById("data-table-BR");
    var tbodyNP = npTable.getElementsByTagName("tbody")[0];
    var tbodyBR = brTable.getElementsByTagName("tbody")[0];

    while (tbodyNP.rows.length > 0) {
      tbodyNP.deleteRow(0);
    }
    while (tbodyBR.rows.length > 0) {
      tbodyBR.deleteRow(0);
    }
}

// SAMPLE INPUT LISTENERS:
document.getElementById('octaFrame').addEventListener('click', function() {
        // First we clear out any data inside the NP and BR tables and run the replot function:
        clearBRnNPTables()
        update3DPlot()

        // Next we assign values to max_x, max_y, and max_z
        document.getElementById('maxX').value = 50;
        document.getElementById('maxY').value = 50;
        document.getElementById('maxZ').value = 50;

        // Then apply inputs to the BR and NP Tables:
        var table = document.getElementById("data-table-BR");
        var tbody = table.getElementsByTagName("tbody")[0];

        // Create a new row
        /// First disable buttons:
        document.getElementById('maxX').disabled = true;
        document.getElementById('maxY').disabled = true;
        document.getElementById('maxZ').disabled = true;
        document.querySelector('.disabled-text').style.display = 'inline'


        // Insert many a new row manually:
        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(25, 25, 0)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = ""; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = ""; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = ""; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(25, 25, 50)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = ""; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = ""; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = ""; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(0, 25, 25)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = ""; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = ""; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = ""; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(50, 25, 25)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = ""; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = ""; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = ""; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(25, 0, 25)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = ""; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = ""; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = ""; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(25, 50, 25)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = ""; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = ""; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = ""; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        // Lastly populate the 3D Plot:
        update3DPlot()

    });

document.getElementById('asymObject').addEventListener('click', function() {
        // First we clear out any data inside the NP and BR tables and run the replot function:
        clearBRnNPTables()
        update3DPlot()

        // Next we assign values to max_x, max_y, and max_z
        document.getElementById('maxX').value = 50;
        document.getElementById('maxY').value = 50;
        document.getElementById('maxZ').value = 50;
        document.getElementById('maxBPs').value = 10000;


        // Then apply inputs to the BR and NP Tables:
        var table = document.getElementById("data-table-BR");
        var tbody = table.getElementsByTagName("tbody")[0];

        // Create a new row
        /// First disable buttons:
        document.getElementById('maxX').disabled = true;
        document.getElementById('maxY').disabled = true;
        document.getElementById('maxZ').disabled = true;
        document.querySelector('.disabled-text').style.display = 'inline'

        // Insert many a new row manually:
        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(25, 0, 25)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = ""; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = ""; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = ""; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(15, 50, 25)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = "(35, 50, 25)"; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = ""; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = ""; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(35, 35, 50)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = "(35, 15, 50)"; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = "(15, 15, 50)"; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = "(15, 35, 50)"; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(25, 25, 0)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = ""; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = ""; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = ""; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var table = document.getElementById("data-table-NP");
        var tbody = table.getElementsByTagName("tbody")[0];
        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(0, 15, 15)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = "(20, 35, 35)"; // Replace with your desired value

        var cell3 = newRow.insertCell(2);
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell3.appendChild(removeButton);

        // Lastly populate the 3D Plot:
        update3DPlot()

    });

document.getElementById('wireframeBeam').addEventListener('click', function() {
        // First we clear out any data inside the NP and BR tables and run the replot function:
        clearBRnNPTables()
        update3DPlot()

        // Next we assign values to max_x, max_y, and max_z
        document.getElementById('maxX').value = 70;
        document.getElementById('maxY').value = 37;
        document.getElementById('maxZ').value = 37;

        // Create a new row
        /// First disable buttons:
        document.getElementById('maxX').disabled = true;
        document.getElementById('maxY').disabled = true;
        document.getElementById('maxZ').disabled = true;
        document.querySelector('.disabled-text').style.display = 'inline'


        // Insert many a new row manually:
        // Then apply inputs to the BR and NP Tables:
        var table = document.getElementById("data-table-BR");
        var tbody = table.getElementsByTagName("tbody")[0];

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(0, 28.5, 28.5)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = "(0, 7, 28.5)"; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = "(0, 7, 7)"; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = "(0, 28.5, 7)"; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        var newRow = tbody.insertRow(tbody.rows.length);
        // Create and add cells with specific values
        var cell1 = newRow.insertCell(0);
        cell1.innerHTML = "(70, 28.5, 28.5)"; // Replace with your desired value
        var cell2 = newRow.insertCell(1);
        cell2.innerHTML = "(70, 7, 28.5)"; // Replace with your desired value
        var cell3 = newRow.insertCell(2);
        cell3.innerHTML = "(70, 7, 7)"; // Replace with your desired value
        var cell4 = newRow.insertCell(3);
        cell4.innerHTML = "(70, 28.5, 7)"; // Replace with your desired value
        var cell5 = newRow.insertCell(4);
        // Add the remove button:
        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'remove-br-button';
        removeButton.innerHTML = 'Remove';
        cell5.appendChild(removeButton);

        update3DPlot()

    });
});



