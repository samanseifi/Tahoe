
<B>Tahoe</B>

For public use:

This is a pre-built version of Tahoe configurated for developing the dielectric elastomer models (DE, DEQ1P0, DEQ1P02D and DEQ1P0Elastocapillary). If you want to start working with Tahoe, first you need to read the README file and then take a look at the sourceforce page: <a href="http://tahoe.sourceforge.net/">http://tahoe.sourceforge.net/</a> in order to set the enviroment. The source files can be downloaded from <a href="http://sourceforge.net/projects/tahoe/">here</a>. After downloading the files for building Tahoe simply execute the following command:

<code>~$: ./tahoe-manager init build</code>

For the first timers preparing Tahoe for a specific application is not as simple as it seems. There are some libraries that need to be installed first. You can read the forum <a href="http://tahoe.sourceforge.net/bb/"> here </a> and learn how to set your own version of Tahoe. The development module needs an access to original repository in cloudeforge.

<B>Pre Built Tahoe (This repository)</B>

This for my own private use:

This is a pre-built version of Tahoe  configurated for developing dielectric elastomer models.

<B>Rebuilding in Linux Ubuntu 14.04</B>

It seems that in Ubuntu 12.04 and 14.04 the best choice for compiler is GNU-GCC4.6. This version of Tahoe here, were built in ubuntu 12.04. You need to rebuild the Tahoe in your machine. First make sure you have GCC 4.6 installed in your unix operating system. If you're using linux with GCC 4.6 it should work. However, if there is another versions of GCC is installed probably more recent versions like, for exeample in ubuntu 14.04 the GCC 4.8 is pre-installed. Building Tahoe with GCC 4.8 is not working properly on linux machine. However, for mac OS X one can download GCC 4.8 or more recent versions from fink and it will work fine. For linux versions GCC 4.6 is needed (for now). You can install GCC 4.6 from "Synaptic Package Manager" or directly from ubuntu repositoroes. Following command:

<code>~$: sudo apt-get install gcc-4.6 gfortran-4.6 cpp-4.6</code>

Now it's neccessary to make some changes in the macro files. Simply go to Tahoe/macros and edit the file GNU-GCC4.6.macros: 

<code>~$: cd Tahoe/macros</code>

<code>~/Tahoe/macros$: gedit GNU-GCC4.6.macros</code>

Now change all gcc, g++ and gfortran to gcc-4.6, g++-4.6 and gfortran-4.6. Now you are able to rebuild Tahoe by going to Tahoe subdirectory and build tahoe: 

<code> ./tahoe-manager init build</code>

<B>Rebuilding in OS X</B>

First you need to install GCC compilers by installing Xcode and the command line extension which it contains the GCC compilers. However, probably you won't be able to build Tahoe (you can try), so it is needed to obtain new sets of GNU compilers from package managements such as Homebrew, MacPorts or fink. So far fink has shown promised in order to compile Tahoe successfully. You only need to change the macro settings to fink gcc in Tahoe configuration file (.tahoe_config)
