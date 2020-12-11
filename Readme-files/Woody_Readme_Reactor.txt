Readme- Reactor
Hello to whoever may be reading this file, please read every section of this file so you understand how the file operates
where its dependencies are, and how to change the parameters.

-First we have our import and model solvers in the first 10 lines, do not change this

-I have then created open lists for each one of the parameters I want to see in the figures that I output,
all of them are labeled with units and what specific values they are looking for inside of the code so it should be fairly
simple to see what is being stored.

-If you want to add sections of what is being stored go to lines 252-260, there are examples of how to get values from the functions
as well as values directly from the document

-Since this is a transient equation we have to input initial conditions for most of the parameters in lines 51-76, you can change these to affect the
code if you like but some variables can be touchy and will mess up the entire run if change to much. All of the values used as the initialization
parameters were found from running at steady state conditions. You can change the height of the core through "Hreflector", the radius through
ID and OD. All units are displayed and there meaning in the code.

-Lines 79-100 is the first function which solves four oridnary differential equations in a two loops, one for each minute, and then on lines 243-248
for each minute cycle. This function solves for the precursor density inside the core, the neutron density, the equilibirum coolant temp, and the equilibrium
fuel temp.

-Lines 102 - 108 is the second function which solves for the outlet temperature of the core and the power using params2. You can find these values on line 252

-lines 115-133 are used to get data for the solar model such as the power produced and the mass flow rate in kg/s. In the excel file the change in the mass flow rate
of the solar power on the cold side of the heat exchanger connected to the TES model will automatically change the mass flow rate of the reactor cold side. These values 
are taken moved into a data frame and then converted into a list.

-Then comes the large loop which performs the amount of iterations for the amount of minutes there are in a day. By my account there is 1440 minutes.

-Then comes if statements from line 180-212, there statements read the excel files that was converted into the list. Dpeneding on the mass flow rate set by the solar
power cold side, the corresponding change in the reactor power will be made by the scaling value n_e. This value works logarithimcally.

-Next come more initialization parameters in lines 214-225, these set the initial values for the reactor core functions.

-Lines 227-233 initialize the normalized values of the neutron population, precursor density, starting fuel temp, and starting coolant temp.

-lines 237-248 then perform the solving of the four ordinary differential equations over a minute period of time and that is chopped into 100 sections.

-Lines 265-274 solve out for the thermo hydraulic processes in the core dependent on the fuel temperature since this a natural circulation system.
Vavg is a solved ordinary differential equation solved for analytically and was implemented in the code whcih is through the loop, from this we are able to derive the mass flow rate through
the core and heat thransfer coefficient we are necessary to find to initialize the next loop for solving the reactor kinetics differential equation model.

-lines 277-279 set the inlet values of the hot side of the heatexchanger equal to the outlet of the core fluids temp, and the inlet and outlet of the cold side.

-lines 284-288 update the thermal expasion coefficient and other thermal properties of the coolant based on the outlet and inlet tempereatures of the core.

-lines 301 - 305 set constant properties for the cold side fluid, this can be change if you want a dependence on the inlet and outlet temps of the fluid in the heat
exchanger

-lines 321-325 set the mass flow rate of the cold side the value grabbed from the excel sheet dependent on the mass flow rate of the solar system, and then
solves for the outlet temperature of the hotside fluid that goes back into the core as the inlet fluid. It also solves for the heat exchanger value of heat transfered.

-lines 341 and 342 do not change, these are the saved csv flies of the mass flow rate that are directly fed into the thermal energy storage model which allows the models to
talk to eachother. Changing this will cutoff communication between the codes.

-The last lines plots all of the values called for in the beggining of the file.


