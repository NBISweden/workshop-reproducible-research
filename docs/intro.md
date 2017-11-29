# Intro

## The biological problem

## The exercises

## For Windows users
Install from https://docs.docker.com/docker-for-windows/install/. Note that for older versions of Windows you should use Docker Toolbox.

Then run
```bash
docker run --rm -it -v c:/my_dir:/home/ scilifelablts/reproducible_research_course_slim
```

Note that the idea is that you should edit files in the mounted  c:/my_dir using an editor in your normal OS, say notepad in windows. The terminal in the container is for running stuff, not editing. (nano is installed as an editor if you really want to, we should probable add vi instead).

### Set up
```bash
git clone https://bitbucket.org/scilifelab-lts/reproducible_research_course.git
cd reproducible_research_course
```

## Take down
How to remove stuff

* ~/.ncbi/ # contains settings file, exists only if vdb-config is run
* ~/ncbi/ # created if caching is not disabled by vdb-config
* ~/.bash_profile # conda line
* ~/.conda/
* ~/.condarc
* ~/.continuum/ ?
* ~/.docker/
