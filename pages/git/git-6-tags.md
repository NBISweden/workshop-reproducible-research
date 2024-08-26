Git allows us to _tag_ commits, _i.e._ give names to specific points in the
history of our project. This can be particularly important for reproducible
research, but also for development projects that want to highlight specific
versions of a software. A tag can be, for example, the version of the repository
that was used for the manuscript submission, the version used during
resubmission, and, most importantly, the version used for the final publication.
The first two examples are mainly useful internally, but the latter is essential
for other researchers to be able to rerun your published analysis.

- Let's assume that the status of the repository as it is now is ready for
  a submission to a journal. It may for example contain the scripts that were
  used to generate the manuscript figures. Let's add a tag:

```bash
git tag "submission1"
```

- We can now list all the tags available in the current repository:

```bash
git tag
```

> **Tip** <br>
> You can use the flag `-a` or `--annotate` to give more detailed information
> about a specific tag, similar to a commit message. This can be quite useful
> when there are many changes that happened, in that it allows you to
> summarise them. You can, for example, do `git tag -a submission1 -m
"Annotation for tag submission1"` to write the annotation along with the
> command (similar to the `-m` flag for committing) or just `git tag -a
submission1` to write the annotation with your default editor. To list all
> your tags along with their annotations you can use _e.g._ `git tag -n10`
> (which will list the first 10 lines of each tag's annotation).

- Let's assume we now got comments from the reviewers, and by fixing
  those we had to update our code. Open `config.yml` and change the line
  `max_reads: 25000` to `max_reads: 50000`. Commit and tag the changes:

```bash
git add config.yml
git commit -m "Increase number of reads"
git tag "revision-1"
```

- Now let's say that the reviewers were happy and the manuscript was
  accepted for publication. Let's immediately add a tag:

```bash
git tag "publication"
```

- A good thing about using tags is that you can easily switch between versions
  of your code. Let's move to the first submission version:

```bash
git switch -d submission1
```

When switching between tags you have to use the `-d` (or `--detach`) flag to the
`switch` command. This is because the `switch` command only switches between
branches as default; switching between tags requires us to be in a _detached_
state. This means that we're not actually on a branch, but rather at a specific
point in the Git history without a branch.

- Open `config.yml` and note that the `max_reads` variable is `25000`! To go
  back to the latest version, run:

```bash
git switch main
```

- Open `config.yml` and see that the value is now `50000`.

> **Tip** <br>
> You can also see the difference between tags in the same way as for
> branches and commits using _e.g._ `git diff <tag1> <tag2>`.

At this point could run `git log --oneline --decorate` to get a condensed
commit history, where you should also be able to see the tagged commits.

> **Quick recap** <br>
> We have now learned how to tag important commits:
>
> - `git tag` adds a tag to a commit.
> - `git switch -d` moves between tags in a similar fashion as between
>   branches.
