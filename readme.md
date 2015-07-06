<B>Tahoe</B>

For public use:

This is a pre-built version of Tahoe specifically configurated for developing the "Dielectric Elastomer" models. So if you want to start with Tahoe you need to start reading README file and also take a look at the sourceforce page: <a href="http://tahoe.sourceforge.net/">http://tahoe.sourceforge.net/</a>. The source files can be downloaded from <a href="http://sourceforge.net/projects/tahoe/">here</a>. After downloading the files for building Tahoe simply execute the following command:

<code>~$: ./tahoe-manager init build</code>

Preparing Tahoe for use is not quite an easy job to do. There are some libraries need to be installed first. For more details on installing modules and other details take a look at README file in this repository.

<B>Pre Built Tahoe</B>

For my own private use:

This is a pre built version of Tahoe specifically configurated for developing Dielectric Elastomer models.

<B>Rebuilding Linux Ubuntu 14.04</B>

It seems for linux (Ubuntu 12.04 and 14.04 is tested) the best choose for compiler is GNU-GCC4.6. This version of Tahoe pre built in ubuntu 12.04. You need to rebuild the Tahoe in your machine. First make sure you have GCC 4.6 installed in your machine. If you're using linux with GCC 4.6 you should be fine. However, if you have other versions of GCC for exeample ubuntu 14.04 uses GCC 4.8 you need to install GCC 4.6 using "Synaptic Package Manager" or directly from ubuntu repository. Following command:

<code>~$: sudo apt-get install gcc-4.6 gfortran-4.6 cpp-4.6</code>

Now it's neccessary to make some changes in the macro files. Simply go to Tahoe/macros and edit the file GNU-GCC4.6.macros: 

<code>~$: cd Tahoe/macros</code>

<code>~/Tahoe/macros$: gedit GNU-GCC4.6.macros</code>

Now change all gcc, g++ and gfortran to gcc-4.6, g++-4.6 and gfortran-4.6. Now you are able to rebuild Tahoe by going to tahoe subdirectory and build tahoe. 

<code> ./tahoe-manager init build</code>
