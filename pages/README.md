# Uploading pages to Canvas

The `upload-page-to-canvas.sh` bash script can be used to upload or update
pages on Canvas. Usage is simple, for example:

```bash
./upload-age-to-canvas.sh conda/conda-1-introduction.md
```

This will either upload or update the `conda-1-introduction.md` page to Canvas,
as applicable, using the default parameters as specified in the script. Things
you may want to change include the `COURSE_ID` parameter and the maximum width
of the resulting page.

## Page names

Canvas uses some automation for page naming by using the original filename.
A page using the `conda-2-the-basics.md` original file will be named `Conda
2 The Basics` at Canvas. You may manually edit this to change capitalisation
or add special characters and subsequent updates to the same page will keep
those changes. For example, manually editing the page in Canvas to be `Conda 2:
The basics` will keep that change after further updates, while changing it to
*e.g.* `Conda 2: Some basics` will cause a new page with the original, unaltered
name to be created.
