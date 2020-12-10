The Steam Generator code is set up as follows

The first section of code (lines 1-3) import varius needed moduels as well as importing the TES code.

The second section of code(lines 5 - 35) is the required inputs of the system. 

	Each variable required has a description of the needed component as well as the units needed.

The third section of code (lines 38 - 122) define each section of the heat exchanger and use the effectiveness - NTU method to calculate the unknown temperatures. The return values of each section are the variables input into the next section as the hot side input temperature. DO NOT CHANGE THIS IT WILL BREAK THE CODE

The fourth section of code (lines 124 - 142) define the effectiveness equation required for the above section. DO NOT CHANGE THIS IT WILL BREAK THE CODE.

The fifth section of code (lines 145 - 152) are the lines that call each definition from the sections of heat exchanger. The required variables are called from the second section. The required inputs for this section are the area as the second to last variable as well as the type of heat exchanger (counter flow, parallel flow, or Double shell pass)

The sixth section of code (lines 159 - 181) are the lines that calcualte the Q values and efficiencys based on the temperatures calcualted from section three, the efficiency of the system based on the Q values of water and salt to make sure they are the same or similar value. Finally the percent energy transfer of each section is calcualted.


