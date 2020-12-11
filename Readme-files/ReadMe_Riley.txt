Riley Willt Python Files

In order to run these files you need to install the IAWPS library for python.
Install iapws (steam tables) via iPython Console or command line <pip install iapws> Or hit the following link: https://pypi.org/project/iapws/

1. ParameterVariation_Riley.py

This file contains load variation, temperature variation, and other parameter variations.
It also runs the whole loop of everyone elses code.
For the load variation, PowerBlock1_Riley.py is used.
For the temperature variation, PowerBlock2_Riley.py is used.
For the other parameter variation, PowerBlock3_Riley.py is used.
Here are links to where I gathered the demand curve data at.
MISO Data - https://www.misoenergy.org/markets-and-operations/real-time--market-data/market-reports/#nt=%2FMarketReportType%3AEIA%2FMarketReportName%3AMISO%20Daily%20Reporting%20(xml)&t=500&p=0&s=MarketReportPublished&sd=desc
CAISO Data - http://www.caiso.com/TodaysOutlook/Pages/default.aspx

# 2. GUI_Riley.py

When this file is ran it displays a GUI of the power block model and the temperature at every state.
It uses PowerBlock_Riley.py to do this.
If you run into an error while running this file, just restart the kernel.

# 3. PowerBlock_Riley.py

Power block cycle For the GUI.

# 4. PowerBlock1_Riley.py

For the load variation in ParameterVariation_Riley.py. Returns different values.

# 5. PowerBlock2_Riley.py

For the temperature variation. Comments out mainTemp and reheatTemp so that they can be varied by the loop.

# 6. PowerBlock3_Riley.py

For the other parameter variation. Returns different values.

# 7. Condenser.py

This file is the condenser file that Chad made. I modified it slightly to run as a function with an input of the T_steam_in that returns the T_steam_out. I use it in many of my files.

# 8. components.py

This file is from the previous group. It contains the functions for all of the components for the power block.

# 9. Power Block.PNG

This is the picture of the power block that is used in the GUI.
