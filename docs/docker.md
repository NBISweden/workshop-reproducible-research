# Docker tutorial

## Introduction

### Quick intro to Docker


## Practical exercise
### Set up
This exercise depends on files from the course BitBucket repo. Take a look at the [intro](index) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/git_jupyter_docker`.

First we need to install Docker. This is quite straightforward on OSX or Windows and a little more cumbersome on Linux. Note that Docker runs as root, which means that you have to have sudo privileges on your computer in order to install or run Docker.

#### OSX
Go to  [https://docs.docker.com/docker-for-mac/install/#download-docker-for-mac](https://docs.docker.com/docker-for-mac/install/#download-docker-for-mac) and select "Get Docker for Mac (Stable)". This will download a dmg file. Click on it once it's done to start the installation. This will open up a window where you can drag the Docker.app to Applications. Close the window and click the Docker app from the Applications menu. Now it's basically just to click "next" a couple of times and we should be good to go. You can find the Docker icon in the menu bar in the upper right part of the screen.

#### Windows
The instructions are different depending on if you have Windows 10 or Windows 7 (earlier versions aren't supported). In order to run Docker on Windows your computer must support Hardware Virtualization Technology and virtualization must be enabled. This is typically done in BIOS. This is outside the scope of this tutorial, so we'll simply go ahead as if though it's enabled and hope that it works.

On Windows 10 we will install Docker for Windows, which is available at [https://docs.docker.com/docker-for-windows/install/#download-docker-for-windows](https://docs.docker.com/docker-for-windows/install/#download-docker-for-windows). Select "Get Docker for Windows (Stable)".

1. Once it's downloaded, double-click Docker for Windows Installer.exe to run the installer.

2. Follow the install wizard and accept the license, authorize the installer, and proceed with the install. You will be asked to authorize Docker.app with your system password during the install process. Click Finish to exit the installer.

3. Start Docker from the Start menu. You can search for it if you cannot find it. The Docker whale icon should appear in the task bar.

4. Now we want to share your local drive(s), so that they are available for Docker. Right-click on the Docker whale icon in the task bar and select "Settings". Go to "Shared drives" and enable the drives you want Docker to have access to. Note that the drive where you'll be running the tutorials from has to be enabled (most likely `C:\`).

On Windows 7 we will instead use Docker Toolbox, which is available at [https://docs.docker.com/toolbox/toolbox_install_windows/](https://docs.docker.com/toolbox/toolbox_install_windows/). Select "Get Docker Toolbox for Windows".

1. Install Docker Toolbox by double-clicking the installer. Step through the installation and accept all the defaults. If Windows security dialog prompts you to allow the program to make a change, choose Yes.

2. You should now have a Docker Quickstart icon on the desktop.


#### Linux
How to install Docker differs a bit depending on your Linux distro, but the steps are the same. For details on how to do it on your distro see [https://docs.docker.com/engine/installation/#server](https://docs.docker.com/engine/installation/#server).

Here we show how to do it for Ubuntu, which is the most common desktop distribution. The same instructions apply to distributions based on Ubuntu, such as Elementary OS or Linux Mint. Docker requires a 64-bit Ubuntu version 14.04 or higher. If your OS is from 2015 or earlier you can double check this with `lsb_release -a`.

1. Add the GPG key for the official Docker repository to the system:
```bash
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
```

2. Add the Docker repository to APT sources:
```bash
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
```

3. Update the package database with the Docker packages from the new repo:
```bash
sudo apt-get update
```

4. Install Docker Community Edition:
```bash
sudo apt-get install -y docker-ce
```

5. Docker should now be installed, the daemon started, and the process enabled to start on boot. Check that it's running:
```bash
sudo systemctl status docker
```
The output should say something about "Active: active (running) since..".

As mentioned before, Docker needs to run as root. You can archive this by prepending all Docker commands with `sudo`. This is the approch that we will take in this tutorial, since the set up becomes a little simpler. If you plan on continuing using Docker you can get rid of this by adding your user to the group `docker`. Here are instructions for how to do this: [https://docs.docker.com/engine/installation/linux/linux-postinstall/](https://docs.docker.com/engine/installation/linux/linux-postinstall/).
