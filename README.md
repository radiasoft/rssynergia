### RadiaSoft Synergia Support
Support tools for [Synergia](https://web.fnal.gov/sites/Synergia).

Read the Docs: http://rssynergia.readthedocs.org/en/latest/

#### License
License: http://www.apache.org/licenses/LICENSE-2.0.html

Copyright (c) 2015 RadiaSoft LLC.  All Rights Reserved.

### Run Synergia on the RadiaSoft JupyterHub server

    Point your browser to the following URL
    https://jupyter.radiasoft.org
    
    Click on the "Sign in with GitHub" button
    A GitHub login is required. If you're not already logged in, then please do so.
    If you have not previously authorized the JupyterHub server on GitHub, please do so now.
    If you see a green "My Server" button, then click on it.
    If you see a Jupyter server home page, then look towards the upper right part of your browser.
    Click the "New" button, then select "Terminal" from the dropdown menu.
    The Jupyter commandline terminal window will open in a new tab.
    
    In the Jupyter terminal window, clone this repo as follows:
    > git clone https://github.com/radiasoft/rssynergia
    
    Go back to the JupyterHub home page (i.e. the URL referenced above)
    You should see that there is now a directory called 'rssynergia'
    Drill down into this directory by clicking on directory names, to:
    .../rssynergia/examples/drift_expansion/
    
    Click on the Jupyter notebook named "sc_drift_expansion.ipynb"
    This opens the notebook.
    Repeatedly hitting "shift-enter" on your keyboard will execute the cells one by one.


### Create a Synergia development environment on the JupyterHub server
    Work in the Jupyter terminal window you opened above.
    You are in a Python2 environment in Fedora 27 linux
    You are username "vagrant" 
    You have sudo privileges with password "vagrant"
    
    By default, you are in directory /home/vagrant/jupyter/
    Move to the 'radiasoft' working directory
    > cd /home/vagrant/src/radiasoft

#### Install and build Synergia2
    Many particle accelerator and radiation codes have been installed.
    You can see the build directories in the local subdirectory '/home/vagrant/src/radiasoft/codes'
    However, these are not suitable for development.
    
    From within the .../radiasoft/ directory, follow the download instructions for the Synergia source code:
    https://cdcvs.fnal.gov/redmine/projects/contract-synergia2/wiki/Download_and_build_the_current_Synergia_release

    Assuming the build completed without errors, cd into the synergia source directory:
    > cd /home/vagrant/src/radiasoft/synergia2-devel/build/synergia2
    
    Invoke the test suite in this directory via:
    > make test
    
    If some of the tests don't pass, then please create an issue here describing the problem, with a screenshot:
    https://github.com/radiasoft/rssynergia/issues

### Create a Synergia development environment on your laptop or desktop

    Install VirtualBox on your computer
    https://www.virtualbox.org/wiki/Downloads

    Install Vagrant on your computer
    https://www.vagrantup.com/docs/installation/

    Create a directory in a linux, unix or cygwin environment
    > mkdir synergia_devel
    > cd synergia_devel

    In this directory, create a text file named "Vagrantfile" with the following contents:

    # -*- mode: ruby -*-
    Vagrant.configure(2) do |config|
      config.vm.box = "radiasoft/beamsim"
      config.vm.hostname = "rs"
      config.ssh.insert_key = false
      config.ssh.forward_x11 = true
      config.vm.network "forwarded_port", guest:8000, host:8000
      config.vm.synced_folder ".", "/vagrant", disabled: false
    end

    Invoke 'vagrant' as follows
    > vagrant up

    Wait several minutes, because the container is large
    ssh into the container as follows
    > vagrant ssh

    You are in a Python2 environment in Fedora 27 linux
    You are username "vagrant" 
    You have sudo privileges with password "vagrant"
    The directory /vagrant is synced with the local directory you created above.
    
    Move to the 'radiasoft' working directory
    > cd src/radiasoft
    
    Now you can follow the "Install and build Synergia2" instructions above
 
### View the Synergia2 source code
    The source code is found in the following subdirectory
    > cd /home/vagrant/src/radiasoft/synergia2-devel/build/synergia2/src/
    
