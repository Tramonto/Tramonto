/* Data structures for parallel RCB

      rcb_dot has 4 required fields as shown below
      other fields can be added by user
      are just carried along as dots migrate to new processors

      examples:

  int       global;                global id # of dot
  int       local;                 local id # (memory loc) of dot before RCB
  int       proc;                  owner of this dot before RCB
*/

			     /* dot to balance on for RCB */
struct rcb_dot
{				/* dot = point in 3-space */
  double x[3];			/* location of dot */
  double weight;		/* weight of dot - if used must be > 0 */
};

struct rcb_tree
{				/* tree of RCB cuts */
  double cut;			/* position of cut */
  int dim;			/* dimension (012) of cut */
};

struct rcb_median
{				/* RCB cut info */
  double lototal, hitotal;	/* weight in each half of active partition */
  double valuelo, valuehi;	/* position of dot(s) nearest to cut */
  double wtlo, wthi;		/* total weight of dot(s) at that position */
  int countlo, counthi;		/* # of dots at that position */
  int proclo, prochi;		/* unique proc who owns a nearest dot */
};

struct rcb_box
{				/* bounding box */
  double lo[3], hi[3];		/* xyz lo/hi bounds */
};
