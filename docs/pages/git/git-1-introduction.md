<iframe id="iframepdf" src="../../../lectures/git/git.pdf" frameborder="0" width="640" height="480" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe> 

Git is a widely used system (both in academia and industry) for version
controlling files and collaborating on code. It is used to track changes in
(text) files, thereby establishing a history of all edits made to each file,
together with short messages about each change and information about who made
it. Git is mainly run from the command line, but there are several tools that
have implemented a graphical user interface to run Git commands.

Using version control for tracking your files, and edits to those, is an
essential step in making your computational research reproducible. A typical Git
workflow consists of:

* Making distinct and related edits to one or several files
* Committing those changes (*i.e.* telling Git to add those edits to the
  history, together with a message about what those changes involve)
* Pushing the commit to a remote repository (*i.e.* syncing your local project
  directory with one in the cloud)

There are many benefits of using Git in your research project:

* You are automatically forced into a more organized way of working, which is
  usually a first step towards reproducibility.
* If you have made some changes to a file and realize that those were probably
  not a good idea after all, it is simple to view exactly what the changes were
  and revert them.
* If there is more than one person involved in the project, Git makes it easy
  to collaborate by tracking all edits made by each person. It will also handle
  any potential conflicting edits.
* Using a cloud-based repository hosting service (the one you push your commits
  to), like *e.g.* [GitHub](https://github.com/) or
  [Bitbucket](https://bitbucket.org/), adds additional features, such as being
  able to discuss the project, comment on edits, or report issues.
* If at some point your project will be published GitHub or Bitbucket (or
  similar) are excellent places to publicly distribute your code. Other
  researchers can then use Git to access the code needed for reproducing your
  results, in exactly the state it was when used for the publication.
* If needed, you can host private repositories on GitHub and Bitbucket as well.
  This may be convenient during an ongoing research project, before it is
  publicly published.

These tutorials will walk you through the basics of using Git as a tool for
reproducible research. The things covered in these tutorials are what you will
be using most of the time in your day-to-day work with Git, but Git has many
more advanced features that might be of use to you.

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](pre-course-setup) for instructions on how to set it up if you haven't
done so already.
