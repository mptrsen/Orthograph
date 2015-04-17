push:
	# update version number
	sed -i -e "s/'[^']\+'/'$(git describe --abbrev=0)'/" Orthograph/Version.pm
	# make commit
	git commit -a
	# update changelog with most recent commit, and amend
	git log --pretty --date=short > ChangeLog
	git commit --all --amend --no-edit
	# push to github
	git push
