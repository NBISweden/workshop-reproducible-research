The following material contains some additional tips and tricks on how to use
Jupyter notebooks. This is not part of the core of the Jupyter material and you
can choose what you want to go through, or skip it entirely.

Here are some useful resources if you want to read more about Jupyter in
general:

* The [Jupyter project site](http://jupyter.org) contains a lot of information
  and inspiration.
* The [Jupyter Notebook documentation](
  https://jupyter-notebook.readthedocs.io/en/stable/).
* A [guide](http://ipywidgets.readthedocs.io/en/stable/index.html) to using
  widgets for creating interactive notebooks.

## Running Jupyter notebooks on a cluster

* Login to Uppmax, making sure to use a specific login node, _e.g._ `rackham1`:

```
ssh <your-user-name>@rackham1.uppmax.uu.se
```

* Create/activate a Conda environment containing `jupyter`, _e.g._:

```
mamba create -n jupyter -c conda-forge jupyter
```

* activate the environment, then run:

```
jupyter notebook --no-browser
```

When the Jupyter server starts up you should see something resembling:
```
[I 2023-11-13 22:15:36.944 ServerApp] Serving notebooks from local directory: <path-to-your-directory>
[I 2023-11-13 22:15:36.944 ServerApp] Jupyter Server 2.10.0 is running at:
[I 2023-11-13 22:15:36.944 ServerApp] http://localhost:8888/tree?token=25fa07e89b7c0bc2e518f259ba79c67847ca813cdf4eeed6
[I 2023-11-13 22:15:36.944 ServerApp]     http://127.0.0.1:8888/tree?token=25fa07e89b7c0bc2e518f259ba79c67847ca813cdf4eeed6
[I 2023-11-13 22:15:36.944 ServerApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
```

Now a Jupyter notebook server is running on the Uppmax end. The line that says:

```
[I 2023-11-13 22:15:36.944 ServerApp] http://localhost:8888/tree?token=25fa07e89b7c0bc2e518f259ba79c67847ca813cdf4eeed6
```

Contains information on the port used on the server side (8888 in this case) and
the token required to use the server (`25fa07e89b7c0bc2e518f259ba79c67847ca813cdf4eeed6`).

Next step is to use this information to login to the server from your local
computer.

**On your local computer**

In a terminal, run the following command to start port forwarding of
port 8080 on your local computer to the remote port on the Uppmax side. Replace
<remote-port> with the port given when you started the server on Uppmax. Also
replace <your-user-name> with your user name on Uppmax.

```bash
ssh -N -L localhost:8080:localhost:<remote-port> <your-user-name>@rackham1.uppmax.uu.se
```

As long as this process is running the port forwarding is running. To disable it
simply interrupt it with `CTRL + C`.

Connect to the Jupyter server by opening `localhost:8080` in your browser. When
prompted, paste the token you got when starting the server on Uppmax and set a
new password.

## Using Binder to share interactive notebooks

[Binder](https://mybinder.org/) is a service that allows you to share Jupyter
notebooks with others, while also allowing them to run the notebooks in the
browser. This is great if you wish to share an analysis and have others interact
with the code and results, without them having to install anything locally. What
you will need is:

1. A public GitHub repository containing the notebooks you want to share.
2. An `environment.yml` file in the repository containing the Conda environment
   required to run the notebooks.
3. Data files (if any) required to run the notebook(s).

Binder will then create a Docker image containing the Conda environment and the
notebooks, and run a Jupyter server on this image. The Docker image is then
hosted on the Binder server and can be used by anyone with the link to the
repository to run the notebooks interactively in their browser.

To show you an example we've created a basic [GitHub
repository](https://github.com/NBISweden/workshop-reproducible-research-binder_example)
containing the `supplementary_material.ipynb` notebook from the previous
section. If you go to the repository you will see a badge saying "launch
binder", click this to start the Binder server. This will take a few minutes the
first time you do it, but after that it should be faster. When the server is
ready you will be presented with the now familiar Jupyter interface. Go ahead
and open up the `supplementary_material.ipynb` notebook and run it.

You can now interact with the notebook as you would if you had it running on a
local Jupyter server. You can change the code, run it, and see the results. You
can also add new cells and write new code. However, you cannot save the changes
you make to the notebook.

To read more about Binder and how to use it, see the 
[Binder documentation](https://mybinder.readthedocs.io/en/latest/). For pointers
on how to make data available to the notebooks you share via Binder, see this
guide on
[Accessing data in your
Binder](https://the-turing-way.netlify.app/communication/binder/zero-to-binder.html#accessing-data-in-your-binder).