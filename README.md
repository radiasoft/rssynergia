### RadiaSoft Synergia Support
Support tools for [Synergia](https://web.fnal.gov/sites/Synergia).

Read the Docs: http://rssynergia.readthedocs.org/en/latest/

#### License
License: http://www.apache.org/licenses/LICENSE-2.0.html

Copyright (c) 2015 RadiaSoft LLC.  All Rights Reserved.

### Create a Synergia development environment

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
    
 
