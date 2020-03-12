#
# Regular cron jobs for the gravity package
#
0 4	* * *	root	[ -x /usr/bin/gravity_maintenance ] && /usr/bin/gravity_maintenance
