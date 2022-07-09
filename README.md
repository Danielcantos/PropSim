# PropSim
Pseudo 1-D performance simulator of solid propellant rocket motors

PropSim simulates the operation of a solid propellant rocket motor by discretizing the propellant grain into finite volumes. It takes into account motor ignition, steady-state and tail-off phases. It considers both ideal and non-ideal phenomena, such as nozzle throat erosion and presence of propellant cracks.

DISCLAIMER: HEAVY WIP, Personal project. Functionality IS NOT guaranteed.

Installation instructions:
- Download the .exe installer. This is an app created via Matlab app designer so this executable will install both the program itself and the Matlab Runtime libraries necessary for its correct operation.
- Install as any other application. A Matlab install window will guide you through the process.
- IMPORTANT: Use the default installation directory for PropSim. That is C:\Program Files\PropSIM. Using any other directory will BREAK the program. A shortcut can be created to execute the program. Work is being done to solve these problems.
- IMPORTANT: The shortcut/app should be executed in administrator mode. It will work without it but if you want to create and save new cases, materials and other important files, write permissions for the C:\Program Files\PropSIM folder are required, which demand administrator mode.
- Output and log will be generated inside of the folder from which the program/shortcut is executed. You can delete them.
