The following extra material contains some more advanced things you can do with
Conda and the command line in general, which is not part of the main course
materials. All the essential skills of Conda are covered by the previous
section: the material here should be considered tips and tricks from people who
use Conda as part of their daily work. You thus don't need to use these things
unless you want to, and you can even skip this part of the lesson if you like!

## Configuring Conda

The behaviour of your Conda installation can be changed using an optional
configuration file `.condarc`. On a fresh Conda install no such file is
included but it's created in your home directory as `~/.condarc` the first time
you run `conda config`.

You can edit the `.condarc` file either using a text editor or by way of the
`conda config` command. To list all config parameters and their settings run:

```bash
conda config --show
```

Similar to Conda environment files, the configuration file is in YAML syntax.
This means that the config file is structured in the form of `key:value` pairs
where the `key` is the name of the config parameter (*e.g.* `auto_update_conda`)
and the `value` is the parameter setting (*e.g.* `True`).

Adding the name of a config parameter to `conda config --show` will show only
that parameter, *e.g.* `conda config --show channels`.

You can change parameters with the `--set`, `--add`, `--append` and `--remove`
flags to `conda config`.

If you for example want to enable the 'Always yes' behaviour which makes Conda
automatically choose the `yes` option, such as when installing, you can run:

```bash
conda config --set always_yes True
```

To see details about a config parameter you can run `conda config --describe
parameter`. Try running it on the `channels` parameter:

```bash
conda config --describe channels
```

In the beginning of this tutorial we added Conda channels to the `.condarc`
file using `conda config --add channels`. To remove one of the channels from
the configuration file you can run:

```bash
conda config --remove channels conda-forge
```

Check your `.condarc` file to see the change. To add the *conda-forge* channel
back to the top of the `channels` simply run:

```bash
conda config --add channels conda-forge
```

To completely remove a parameter and all its values run:

```bash
conda config --remove-key parameter
```

For a list of Conda configuration parameters see the
[Conda configuration](https://docs.conda.io/projects/conda/en/latest/configuration.html)
page.

## Managing Python versions

With Conda it's possible to keep several different versions of Python on your
computer at the same time, and switching between these versions is very easy.
However, a single Conda environment can only contain one version of Python.

### Your current Python installation

The Conda base environment has its own version of Python installed.
When you open a terminal (after having installed Conda on your system) this base
environment is activated by default (as evidenced by `(base)` prepended to your
prompt). You can check what Python version is installed in this environment by
running `python --version`. To see the exact path to the Python executable type
`which python`.

In addition to this your computer may already have Python installed in a
separate (system-wide) location outside of the Conda installation. To see if
that is the case type `conda deactivate` until your prompt is not prepended
with a Conda environment name. Then type `which python`. If a path was printed
to the terminal (*e.g.* `/usr/bin/python`) that means some Python version is
already installed in that location. Check what version it is by typing `python
--version`.

Now activate the base Conda environment again by typing `conda activate` (or
the equivalent `conda activate base`) then check the Python installation path
and version using `which` and `python --version` as above. See the difference?
When you activate a Conda environment your `$PATH` variable is updated so that
when you call `python` (or any other program) the system first searches the
directory of the currently active environment.

### Different Python versions

When you create a new Conda environment you can choose to install a specific
version of Python in that environment as well. As an example, create an
environment containing Python version `3.5` by running:

```bash
conda create -n py35 python=3.5
```

Here we name the environment `py35` but you can choose whatever name you want.

To activate the environment run:

```bash
conda activate py35
```

You now have a completely separate environment with its own Python version.

Let's say you instead want an environment with Python version `2.7` installed.
You may for instance want to run scripts or packages that were written for
Python 2.x and are thus incompatible with Python 3.x. Simply create the new
Conda environment with:

```bash
conda create -n py27 python=2.7
```

Activate this environment with:

```bash
conda activate py27
```

Now, switching between Python versions is as easy as typing `conda activate
py35` / `conda activate py27`.

!!! Note
    If you create an environment where none of the packages require Python,
    *and* you don't explicitly install the `python` package then that new
    environment will use the Python version installed in your base Conda
    environment.

## Decorating your prompt

By default, Conda adds the name of the currently activated environment to the
end of your command line prompt. This is a good thing, as it makes it easier to
keep track of what environment and packages you have access to. The way this is
done in the default implementation becomes an issue when using absolute paths
for environments (specifying `conda env create -p path/to/environment`,
though, as the entire path will be added to the prompt. This can take up a lot
of unnecessary space, but can be solved in a number of ways.

The most straightforward way to solve this is to change the Conda configuration
file, specifically the settings of the `env_prompt` configuration value which
determines how Conda modifies your command line prompt. For more information
about this setting you can run `conda config --describe env_prompt` and to see
your current setting you can run `conda config --show env_prompt`.

By default `env_prompt` is set to `({default_env})` which modifies your prompt
with the active environment name if it was installed using the `-n` flag or if
the environment folder has a parent folder named `envs/`. Otherwise the full
environment path (*i.e.* the 'prefix') is displayed.

If you instead set env_prompt to `({name}) ` Conda will modify your prompt with
the folder name of the active environment. You can change the setting by
running `conda config --set env_prompt '({name}) '`

If you wish to keep the `({default_env})` behaviour, or just don't want to
change your Conda config, an alternative is to keep Conda environment folders
within a parent folder called `envs/`. This will make Conda only add the folder
name of the Conda environment to your prompt when you activate it.

As an example, say you have a project called *project_a* with the project path
`~/myprojects/project_a`. You could then install the environment for *project_a*
into a folder `~/myprojects/project_a/envs/project_a_environment`. Activating
the environment by pointing Conda to it (*e.g.*
`conda activate ~/myprojects/project_a/envs/project_a_environment`) will only
cause your prompt to be modified with *project_a_environment*.

## Bash aliases for conda

Some programmers like to have aliases (_i.e._ shortcuts) for common commands.
Two aliases that might be usefol for you are `alias coac='conda activate'` and
`alias code='conda deactivate'`. Don't forget to add them to your
`~/.bash_profile` if you want to use them!

## Rolling back to an earlier version of the environment

Conda keeps a history of the changes to an environment. You can see
revisions to an environment by using:

```bash
conda list --revisions
```

which shows each revision (numbered) and what's installed.

You can revert back to particular revision using:

```bash
conda install --revision 5
```
