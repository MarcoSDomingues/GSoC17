# GSoC17
My work for Google Summer of Code 2017 (BRL-CAD)

## GPU Boolean Evaluation for CSG Ray-Tracing
During the summer of 2017, I contributed to the BRL-CAD Project (https://brlcad.org/) under Google Summer of Code program. GSoC is a summer program hosted by Google and it is focused in introducing students to open source software development.

My work for BRL-CAD consisted in implementing boolean evaluation for ray-tracing of Constructive Solid Geometry (CSG) in OpenCL.

**Project Plan:**       https://brlcad.org/wiki/User:Marco-domingues/GSoC17/Project              
**Development Logs:**   https://brlcad.org/wiki/User:Marco-domingues/GSoC17/Log


Here you can see the final version of the changes I made during the program. 
The code was applied to BRL-CAD repository: https://sourceforge.net/projects/brlcad/

___
```diff

Index: include/rt/region.h
===================================================================
--- include/rt/region.h	(revision 70095)
+++ include/rt/region.h	(working copy)
@@ -65,6 +65,40 @@
 #define REGION_NULL     ((struct region *)0)
 #define RT_CK_REGION(_p) BU_CKMAG(_p, RT_REGION_MAGIC, "struct region")
 
+#ifdef USE_OPENCL
+struct cl_bool_region {
+    cl_uint btree_offset;           /**< @brief index to the start of the bit tree */
+    cl_int reg_aircode;             /**< @brief Region ID AIR code */
+    cl_int reg_bit;                 /**< @brief constant index into Regions[] */
+    cl_short reg_all_unions;        /**< @brief 1=boolean tree is all unions */
+};
+
+/*
+ * Values for Shader Function ID.
+ */
+#define SH_NONE 	0
+#define SH_PHONG	1
+
+
+struct cl_phong_specific {
+    cl_double wgt_specular;
+    cl_double wgt_diffuse;
+    cl_int shine;
+};
+
+/**
+ * The region structure.
+ */
+struct cl_region {
+    cl_float color[3];		/**< @brief explicit color:  0..1  */
+    cl_int mf_id;
+    union {
+	struct cl_phong_specific phg_spec;
+    }udata;
+};
+
+#endif
+
 /* Print a region */
 RT_EXPORT extern void rt_pr_region(const struct region *rp);
 
Index: include/rt/rt_instance.h
===================================================================
--- include/rt/rt_instance.h	(revision 70095)
+++ include/rt/rt_instance.h	(working copy)
@@ -359,6 +359,12 @@
 RT_EXPORT extern void
 clt_db_store_bvh(size_t count, struct clt_linear_bvh_node *nodes);
 
+RT_EXPORT extern void
+clt_db_store_regions(size_t sz_btree_array, struct bit_tree *btp, size_t nregions, struct cl_bool_region *regions, struct cl_region *mtls);
+
+RT_EXPORT extern void
+clt_db_store_regions_table(cl_uint *regions_table, size_t regions_table_size);
+
 RT_EXPORT extern void clt_db_release(void);
 
 
Index: include/rt/shoot.h
===================================================================
--- include/rt/shoot.h	(revision 70095)
+++ include/rt/shoot.h	(working copy)
@@ -176,12 +176,24 @@
     cl_uint seg_sti;
 };
 
+struct cl_partition {
+    struct cl_hit inhit;
+    struct cl_hit outhit;
+    cl_uint inseg;
+    cl_uint outseg;
+    cl_uint forw_pp;                /* index to the next partition */
+    cl_uint back_pp;                /* index to the previous partition */
+    cl_uint region_id;              /* id of the "owning" region */
+    cl_char inflip;                 /* flip inhit->hit_normal */
+    cl_char outflip;                /* flip outhit->hit_normal */
+};
+
 RT_EXPORT extern void
 clt_frame(void *pixels, uint8_t o[3], int cur_pixel, int last_pixel,
 	  int width, int ibackground[3], int inonbackground[3],
 	  double airdensity, double haze[3], fastf_t gamma,
           mat_t view2model, fastf_t cell_width, fastf_t cell_height,
-          fastf_t aspect, int lightmodel);
+          fastf_t aspect, int lightmodel, int a_no_booleans);
 #endif
 
 
Index: include/rt/tree.h
===================================================================
--- include/rt/tree.h	(revision 70095)
+++ include/rt/tree.h	(working copy)
@@ -254,6 +254,49 @@
 
 #define TREE_LIST_NULL  ((struct tree_list *)0)
 
+#ifdef USE_OPENCL
+/**
+ * Flattened version of the infix union tree.
+ */
+#define UOP_UNION        1         /**< @brief  Binary: L union R */
+#define UOP_INTERSECT    2         /**< @brief  Binary: L intersect R */
+#define UOP_SUBTRACT     3         /**< @brief  Binary: L subtract R */
+#define UOP_XOR          4         /**< @brief  Binary: L xor R, not both*/
+#define UOP_NOT          5         /**< @brief  Unary:  not L */
+#define UOP_GUARD        6         /**< @brief  Unary:  not L, or else! */
+#define UOP_XNOP         7         /**< @brief  Unary:  L, mark region */
+
+#define UOP_SOLID        0         /**< @brief  Leaf:  tr_stp -> solid */
+
+/**
+ * bit expr tree representation
+ *
+ * node:
+ *      uint uop : 3
+ *      uint right_child : 29
+ *
+ * leaf:
+ *      uint uop : 3
+ *      uint st_bit : 29
+ */
+struct bit_tree {
+    uint val;
+};
+
+struct cl_tree_bit {
+    cl_uint val;
+};
+
+/* Print a bit expr tree */
+RT_EXPORT extern void rt_pr_bit_tree(const struct bit_tree *btp,
+                                     int idx,
+                                     int lvl);
+
+RT_EXPORT extern void rt_bit_tree(struct bit_tree *btp,
+                                  const union tree *tp,
+                                  size_t *len);
+#endif
+
 /* Print an expr tree */
 RT_EXPORT extern void rt_pr_tree(const union tree *tp,
 	                         int lvl);
Index: src/librt/CMakeLists.txt
===================================================================
--- src/librt/CMakeLists.txt	(revision 70095)
+++ src/librt/CMakeLists.txt	(working copy)
@@ -208,7 +208,7 @@
   primitives/superell/superell.c
   primitives/superell/superell_brep.cpp
   primitives/superell/superell_mirror.c
-  primitives/table.c
+  primitives/table.cpp
   primitives/tgc/tgc.c
   primitives/tgc/tgc_brep.cpp
   primitives/tgc/tgc_mirror.c
@@ -276,6 +276,7 @@
   primitives/revolve/revolve.h
   primitives/rt.cl
   primitives/solver.cl
+  primitives/bool.cl
   primitives/sph/benchmark.sh
   primitives/sph/sph_shot.cl
   primitives/tgc/tgc_shot.cl
@@ -300,6 +301,7 @@
   primitives/common.cl
   primitives/rt.cl
   primitives/solver.cl
+  primitives/bool.cl
 
   primitives/arb8/arb8_shot.cl
   primitives/bot/bot_shot.cl
Index: src/librt/pr.c
===================================================================
--- src/librt/pr.c	(revision 70095)
+++ src/librt/pr.c	(working copy)
@@ -624,7 +624,67 @@
     if (lvl == 0) bu_log("\n");
 }
 
+#ifdef USE_OPENCL
+/**
+ * Produce representations of this bit bool tree
+ */
+void
+rt_pr_bit_tree(const struct bit_tree *btp, int idx, int lvl)
+/* Tree to print */
+/* Offset in tree */
+/* Recursion level */
+{
+    uint uop, val;
 
+    uop = btp[idx].val & 7;
+    val = btp[idx].val >> 3;
+
+    if (lvl == 0) bu_log("bit tree: ");
+
+    switch (uop) {
+        case UOP_SOLID:
+            /* Tree leaf */
+            bu_log("%ld", val);
+            if (lvl == 0) bu_log("\n");
+            return;
+        case UOP_SUBTRACT:
+            bu_log("(");
+            rt_pr_bit_tree(btp, idx+1, lvl+1);
+            bu_log(" %c ", DB_OP_SUBTRACT);
+            rt_pr_bit_tree(btp, val, lvl+1);
+            bu_log(")");
+            break;
+        case UOP_UNION:
+            bu_log("(");
+            rt_pr_bit_tree(btp, idx+1, lvl+1);
+            bu_log(" %c ", DB_OP_UNION);
+            rt_pr_bit_tree(btp, val, lvl+1);
+            bu_log(")");
+            break;
+        case UOP_INTERSECT:
+            bu_log("(");
+            rt_pr_bit_tree(btp, idx+1, lvl+1);
+            bu_log(" %c ", DB_OP_INTERSECT);
+            rt_pr_bit_tree(btp, val, lvl+1);
+            bu_log(")");
+            break;
+        case UOP_XOR:
+            bu_log("(");
+            rt_pr_bit_tree(btp, idx+1, lvl+1);
+            bu_log(" XOR ");
+            rt_pr_bit_tree(btp, val, lvl+1);
+            bu_log(")");
+            break;
+
+        default:
+            bu_log("rt_pr_bit_tree: bad op[%d]\n", uop);
+            exit(1);
+            break;
+    }
+    if (lvl == 0) bu_log("\n");
+}
+#endif
+
 void
 rt_pr_fallback_angle(struct bu_vls *str, const char *prefix, const double *angles)
 {
Index: src/librt/prep.c
===================================================================
--- src/librt/prep.c	(revision 70095)
+++ src/librt/prep.c	(working copy)
@@ -41,7 +41,10 @@
 #include "raytrace.h"
 #include "bn/plot3.h"
 
+#include "optical.h"
+#include "optical/plastic.h"
 
+
 extern void rt_ck(struct rt_i *rtip);
 
 HIDDEN void rt_solid_bitfinder(register union tree *treep, struct region *regp, struct resource *resp);
@@ -445,6 +448,48 @@
 
 
 #ifdef USE_OPENCL
+static void
+rt_btree_translate(struct rt_i *rtip, struct soltab **primitives, struct bit_tree *btp, size_t start, size_t end, const long n_primitives)
+{
+    size_t i;
+    long j;
+    uint uop, st_bit;
+
+    RT_CK_RTI(rtip);
+
+    for (i=start; i<end; i++) {
+        uop = btp[i].val & 7;
+	if (uop == UOP_SOLID) {
+	    st_bit = btp[i].val >> 3;
+	    for (j = 0; j < n_primitives; j++) {
+		if (st_bit == primitives[j]->st_bit) {
+		    btp[i].val = (rtip->rti_Solids[j]->st_bit << 3) | UOP_SOLID;
+		    break;
+		}
+	    }
+	}
+    }
+}
+
+static void
+build_regions_table(cl_uint *regions_table, struct bit_tree *btp, size_t start, size_t end, const long n_primitives, const size_t n_regions, const long reg_id)
+{
+    size_t i;
+    uint uop, st_bit;
+    uint rt_index;
+
+    rt_index = n_regions/32 + 1;
+    for (i=start; i<end; i++) {
+        uop = btp[i].val & 7;
+	if (uop == UOP_SOLID) {
+            st_bit = btp[i].val >> 3;
+	    if (st_bit < n_primitives) {
+	        regions_table[st_bit * rt_index + (reg_id >> 5)] |= 1 << (reg_id & 31);
+	    }
+        }
+    }
+}
+
 void
 clt_prep(struct rt_i *rtip)
 {
@@ -453,6 +498,7 @@
 
     struct soltab **primitives;
     long n_primitives;
+    size_t n_regions;
 
     RT_CK_RTI(rtip);
 
@@ -518,6 +564,101 @@
 	}
 
 	clt_db_store(n_primitives, primitives);
+
+	n_regions = rtip->nregions;
+
+	if (n_regions != 0) {
+	    /* Build boolean regions */
+	    struct region *regp;
+	    struct cl_bool_region *regions;
+	    struct cl_region *mtls;
+	    struct bit_tree *btree;
+	    cl_uint *regions_table;
+	    size_t sz_regions_table;
+	    size_t sz_btree_array;
+	    size_t len;
+
+	    regions = (struct cl_bool_region*)bu_calloc(n_regions, sizeof(*regions), "regions");
+	    mtls = (struct cl_region*)bu_calloc(n_regions, sizeof(*mtls), "mtls");
+
+	    /* Determine the size of all trees to build one array containing
+	     * the bit trees from all regions.
+	     */
+	    sz_btree_array = 0;
+
+	    i = 0;
+	    for (BU_LIST_FOR(regp, region, &(rtip->HeadRegion))) {
+		const struct mfuncs *mfp;
+		const cl_float unset[3] = {1.0f, 1.0f, 1.0f};
+
+		RT_CK_REGION(regp);
+
+		len = 0;
+		rt_bit_tree(NULL, regp->reg_treetop, &len);
+		sz_btree_array += len;
+
+
+		VMOVE(mtls[i].color, unset);
+		mtls[i].mf_id = SH_PHONG;
+
+		if (regp->reg_mater.ma_color_valid) {
+		    VMOVE(mtls[i].color, regp->reg_mater.ma_color);
+		    mtls[i].mf_id = SH_PHONG;
+		}
+
+		mfp = (const struct mfuncs*)regp->reg_mfuncs;
+		if (mfp) {
+		    if (bu_strcmp(mfp->mf_name, "default") ||
+			    bu_strcmp(mfp->mf_name, "phong") ||
+			    bu_strcmp(mfp->mf_name, "plastic") ||
+			    bu_strcmp(mfp->mf_name, "mirror") ||
+			    bu_strcmp(mfp->mf_name, "glass")) {
+			struct phong_specific *src =
+			    (struct phong_specific*)regp->reg_udata;
+			struct cl_phong_specific *dst =
+			    &mtls[i].udata.phg_spec;
+
+			dst->shine = src->shine;
+			dst->wgt_diffuse = src->wgt_diffuse;
+			dst->wgt_specular = src->wgt_specular;
+
+			mtls[i].mf_id = SH_PHONG;
+		    } else {
+			bu_log("Unknown OCL shader: %s\n", mfp->mf_name);
+		    }
+		}
+
+		i++;
+	    }
+
+	    sz_regions_table = n_primitives * ((n_regions/32) + 1);
+	    btree = (struct bit_tree *)bu_calloc(sz_btree_array, sizeof(struct bit_tree), "region btree array");
+	    regions_table = (cl_uint*)bu_calloc(sz_regions_table, sizeof(cl_uint), "regions_table");
+
+	    len = 0;
+	    i = 0;
+	    for (BU_LIST_FOR(regp, region, &(rtip->HeadRegion))) {
+		RT_CK_REGION(regp);
+
+                regions[i].btree_offset = len;
+		regions[i].reg_aircode = regp->reg_aircode;
+		regions[i].reg_bit = regp->reg_bit;
+		regions[i].reg_all_unions = regp->reg_all_unions;
+
+		rt_bit_tree(btree, regp->reg_treetop, &len);
+		rt_btree_translate(rtip, primitives, btree, regions[i].btree_offset, len, n_primitives);
+		build_regions_table(regions_table, btree, regions[i].btree_offset, len, n_primitives, n_regions, i);
+		i++;
+	    }
+
+	    clt_db_store_regions(sz_btree_array, btree, n_regions, regions, mtls);
+	    clt_db_store_regions_table(regions_table, sz_regions_table);
+	    bu_free(mtls, "mtls");
+	    bu_free(regions, "regions");
+	    bu_free(btree, "region btree array");
+	    bu_free(regions_table, "regions_table");
+	}
+
 	bu_free(primitives, "ordered primitives");
     }
 }
Index: src/librt/primitives/bool.cl
===================================================================
--- src/librt/primitives/bool.cl	(nonexistent)
+++ src/librt/primitives/bool.cl	(working copy)
@@ -0,0 +1,1027 @@
+#include "common.cl"
+
+#if !RT_SINGLE_HIT
+
+#define BOOL_STACKSIZE	128
+
+/**
+ * Flattened version of the infix union tree.
+ */
+#define UOP_UNION        1         /* Binary: L union R */
+#define UOP_INTERSECT    2         /* Binary: L intersect R */
+#define UOP_SUBTRACT     3         /* Binary: L subtract R */
+#define UOP_XOR          4         /* Binary: L xor R, not both*/
+#define UOP_NOT          5         /* Unary:  not L */
+#define UOP_GUARD        6         /* Unary:  not L, or else! */
+#define UOP_XNOP         7         /* Unary:  L, mark region */
+
+#define UOP_SOLID        0         /* Leaf:  tr_stp -> solid */
+
+/* Boolean values.  Not easy to change, but defined symbolically */
+#define BOOL_FALSE	0
+#define BOOL_TRUE	1
+
+/* Operations on dynamic bitarrays */
+inline uint
+bindex(const uint b)
+{
+    return (b >> 5);
+}
+
+inline uint
+bmask(const uint b)
+{
+    return (1 << (b & 31));
+}
+
+inline uint
+isset(global uint *bitset, const uint offset, const uint b)
+{
+    return (bitset[offset + bindex(b)] & bmask(b));
+}
+
+inline uint
+clr(global uint *bitset, const uint offset, const uint b)
+{
+    return (bitset[offset + bindex(b)] &= ~bmask(b));
+}
+
+inline uint
+set(global uint *bitset, const uint offset, const uint b)
+{
+    return (bitset[offset + bindex(b)] |= bmask(b));
+}
+
+/**
+ * When duplicating partitions, the bitarray representing the segments
+ * of the partition also has to be copied
+ */
+inline void
+copy_bv(global uint *bitset, const uint bv_index, const uint copy_to, const uint copy_from)
+{
+    for (uint i = 0; i < bv_index; i++) {
+        bitset[copy_to + i] = bitset[copy_from + i];
+    }
+}
+
+/**
+ * Update 'back_pp' and 'forw_pp' values when inserting new partitions
+ * Update head value when inserting at head
+ *
+ * Head partition: 'back_pp' = head
+ * Tail partition: 'forw_pp' = UINT_MAX
+ */
+inline void
+insert_partition_pp(global struct partition *partitions, int pp_count, uint *head, uint new, uint old)
+{
+    if (pp_count == 0)
+        return;
+
+    if (*head == old) {
+        partitions[old].back_pp = new;
+        partitions[new].back_pp = new;
+        partitions[new].forw_pp = old;
+        *head = new;
+    } else {
+	partitions[partitions[old].back_pp].forw_pp = new;
+	partitions[new].back_pp = partitions[old].back_pp;
+	partitions[new].forw_pp = old;
+	partitions[old].back_pp = new;
+    }
+}
+
+/**
+ * Update 'back_pp' and 'forw_pp' values when appending new partitions
+ * Update tail value
+ *
+ * Head partition: 'back_pp' = head
+ * Tail partition: 'forw_pp' = UINT_MAX
+ */
+inline void
+append_partition_pp(global struct partition *partitions, int pp_count, uint new, uint *tail)
+{
+    if (pp_count == 0) {
+        partitions[new].back_pp = new;
+        partitions[new].forw_pp = UINT_MAX;
+        *tail = new;
+    } else {
+	partitions[new].back_pp = *tail;
+	partitions[new].forw_pp = UINT_MAX;
+	partitions[*tail].forw_pp = new;
+	*tail = new;
+    }
+}
+
+inline void
+initialize_partition(global struct partition *partitions, const uint pp_idx)
+{
+    if (pp_idx != UINT_MAX) {
+        partitions[pp_idx].inflip = 0;
+        partitions[pp_idx].outflip = 0;
+        partitions[pp_idx].region_id = UINT_MAX;
+    }
+}
+
+/**
+ * If a zero thickness segment abuts another partition, it will be
+ * fused in, later.
+ *
+ * If it is free standing, then it will remain as a zero thickness
+ * partition, which probably signals going through some solid an odd
+ * number of times, or hitting an NMG wire edge or NMG lone vertex.
+ */
+void
+bool_weave0seg(RESULT_TYPE segp, global struct partition *partitions, int pp_count, global uint *h, global uint *segs_bv, const uint bv_index, uint k, size_t id, uint start_index, uint *head)
+{
+    global struct partition *pp;
+    global struct partition *newpp;
+
+    //bool_weave0seg() with empty partition list
+    if (pp_count == 0)
+	return;
+
+    /* See if this segment ends before start of first partition */
+    if (segp->seg_out.hit_dist < partitions[*head].inhit.hit_dist) {
+	newpp = &partitions[start_index + pp_count];
+        initialize_partition(partitions, start_index + pp_count);
+
+	newpp->inseg = k;
+	newpp->inhit = segp->seg_in;
+	newpp->outseg = k;
+	newpp->outhit = segp->seg_out;
+	set(segs_bv, (start_index + pp_count) * bv_index, k-h[id]);
+	insert_partition_pp(partitions, pp_count, head, start_index + pp_count, *head);
+	pp_count++;
+	return;
+    }
+
+    /*
+     * Cases: seg at start of pt, in middle of pt, at end of pt, or
+     * past end of pt but before start of next pt.
+     *
+     * XXX For the first 3 cases, we might want to make a new 0-len pt,
+     * XXX especially as the NMG ray-tracer starts reporting wire hits.
+     */
+    for (uint current_index = *head; current_index != UINT_MAX; current_index = partitions[current_index].forw_pp) {
+	pp = &partitions[current_index];
+	if (NEAR_EQUAL(segp->seg_in.hit_dist, pp->inhit.hit_dist, rti_tol_dist) ||
+	    NEAR_EQUAL(segp->seg_out.hit_dist, pp->inhit.hit_dist, rti_tol_dist)
+	   )
+	    return;
+
+	if (NEAR_EQUAL(segp->seg_in.hit_dist, pp->outhit.hit_dist, rti_tol_dist) ||
+	    NEAR_EQUAL(segp->seg_out.hit_dist, pp->outhit.hit_dist, rti_tol_dist)
+	   )
+	    return;
+
+	if (segp->seg_out.hit_dist <= pp->outhit.hit_dist &&
+	    segp->seg_in.hit_dist >= pp->inhit.hit_dist)
+	    return;
+
+	if (pp->forw_pp != UINT_MAX && segp->seg_out.hit_dist < partitions[pp->forw_pp].inhit.hit_dist) {
+	    //0-len segment after existing partition, but before next partition.
+	    newpp = &partitions[start_index + pp_count];
+            initialize_partition(partitions, start_index + pp_count);
+
+	    newpp->inseg = k;
+	    newpp->inhit = segp->seg_in;
+	    newpp->outseg = k;
+	    newpp->outhit = segp->seg_out;
+	    set(segs_bv, (start_index + pp_count) * bv_index, k-h[id]);
+	    insert_partition_pp(partitions, pp_count, head, start_index + pp_count, pp->forw_pp);
+	    pp_count++;
+	    return;
+	}
+    }
+}
+
+__kernel void
+rt_boolweave(global struct partition *partitions, global uint *head_partition, RESULT_TYPE segs,
+        global uint *h, global uint *segs_bv, const int cur_pixel,
+        const int last_pixel, const int max_depth)
+{
+    const size_t id = get_global_size(0)*get_global_id(1)+get_global_id(0);
+
+    if (id >= (last_pixel-cur_pixel+1))
+	return;
+
+    const int pixelnum = cur_pixel+id;
+
+    global struct partition *pp;
+    double diff, diff_se;
+
+    head_partition[id] = UINT_MAX;
+
+    uint start_index = 2 * h[id];
+    uint head_pp = start_index;
+    uint tail_pp = start_index;
+    uint bv_index = max_depth/32 + 1;
+    int pp_count;
+
+    pp_count = 0;
+    for (uint k=h[id]; k!=h[id+1]; k++) {
+	RESULT_TYPE segp = segs+k;
+
+	global struct partition *newpp;
+	uint lastseg;
+	global struct hit *lasthit;
+	bool lastflip = 0;
+	uint j;
+
+	/* Make nearly zero be exactly zero */
+	if (NEAR_ZERO(segp->seg_in.hit_dist, rti_tol_dist))
+	    segp->seg_in.hit_dist = 0;
+	if (NEAR_ZERO(segp->seg_out.hit_dist, rti_tol_dist))
+	    segp->seg_out.hit_dist = 0;
+
+	/* Totally ignore things behind the start position */
+	if (segp->seg_in.hit_dist < -10.0)
+	    continue;
+
+	if (!(segp->seg_in.hit_dist >= -INFINITY && segp->seg_out.hit_dist <= INFINITY))
+	    continue;
+
+	if (segp->seg_in.hit_dist > segp->seg_out.hit_dist)
+	    continue;
+
+	diff = segp->seg_in.hit_dist - segp->seg_out.hit_dist;
+
+	/*
+	 * Weave this segment into the existing partitions, creating
+	 * new partitions as necessary.
+	 */
+	if (pp_count == 0) {
+	    /* No partitions yet, simple! */
+	    pp = &partitions[start_index + pp_count];
+            initialize_partition(partitions, start_index + pp_count);
+
+	    pp->inseg = k;
+	    pp->inhit = segp->seg_in;
+	    pp->outseg = k;
+	    pp->outhit = segp->seg_out;
+	    set(segs_bv, (start_index + pp_count) * bv_index, k-h[id]);
+	    append_partition_pp(partitions, pp_count, start_index + pp_count, &tail_pp);
+	    pp_count++;
+	} else if (NEAR_ZERO(diff, rti_tol_dist)) {
+	    /* Check for zero-thickness segment, within tol */
+	    bool_weave0seg(segp, partitions, pp_count, h, segs_bv, bv_index, k, id, start_index, &head_pp);
+	} else if (pp_count > 0 && segp->seg_in.hit_dist >= partitions[tail_pp].outhit.hit_dist) {
+	    /*
+	     * Segment starts exactly at last partition's end, or
+	     * beyond last partitions end.  Make new partition.
+	     */
+	    pp = &partitions[start_index + pp_count];
+            initialize_partition(partitions, start_index + pp_count);
+
+	    pp->inseg = k;
+	    pp->inhit = segp->seg_in;
+	    pp->outseg = k;
+	    pp->outhit = segp->seg_out;
+	    set(segs_bv, (start_index + pp_count) * bv_index, k-h[id]);
+	    append_partition_pp(partitions, pp_count, start_index + pp_count, &tail_pp);
+	    pp_count++;
+	} else {
+	    /* Loop through current partition list weaving the current
+	     * input segment into the list. The following three variables
+	     * keep track of the current starting point of the input
+	     * segment. The starting point of the segment moves to higher
+	     * hit_dist values (as it is woven in) until it is entirely
+	     * consumed.
+	     */
+	    lastseg = k;
+	    lasthit = &segp->seg_in;
+	    lastflip = 0;
+	    for (j = head_pp; j != UINT_MAX; j = partitions[j].forw_pp) {
+		pp = &partitions[j];
+		diff_se = lasthit->hit_dist - pp->outhit.hit_dist;
+
+		if (diff_se > rti_tol_dist) {
+		    /* Seg starts beyond the END of the
+		     * current partition.
+		     *	PPPP
+		     *	    SSSS
+		     * Advance to next partition.
+		     */
+		    continue;
+		}
+		diff = lasthit->hit_dist - pp->inhit.hit_dist;
+		if (diff_se > -(rti_tol_dist) && diff > rti_tol_dist) {
+		    /*
+		     * Seg starts almost "precisely" at the
+		     * end of the current partition.
+		     *	PPPP
+		     *	    SSSS
+		     * FUSE an exact match of the endpoints,
+		     * advance to next partition.
+		     */
+		    lasthit->hit_dist = pp->outhit.hit_dist;
+		    continue;
+		}
+
+		/*
+		 * diff < ~~0
+		 * Seg starts before current partition ends
+		 *	PPPPPPPPPPP
+		 *	  SSSS...
+		 */
+		if (diff >= rti_tol_dist) {
+		    /*
+		     * lasthit->hit_dist > pp->pt_inhit->hit_dist
+		     * pp->pt_inhit->hit_dist < lasthit->hit_dist
+		     *
+		     * Segment starts after partition starts,
+		     * but before the end of the partition.
+		     * Note:  pt_seglist will be updated in equal_start.
+		     *	PPPPPPPPPPPP
+		     *	     SSSS...
+		     *	newpp|pp
+		     */
+		    /* new partition is the span before seg joins partition */
+		    newpp = &partitions[start_index + pp_count];
+		    *newpp = *pp;
+		    copy_bv(segs_bv, bv_index, (start_index + pp_count) * bv_index, j * bv_index);
+
+		    pp->inseg = k;
+		    pp->inhit = segp->seg_in;
+		    pp->inflip = 0;
+		    newpp->outseg = k;
+		    newpp->outhit = segp->seg_in;
+		    newpp->outflip = 1;
+		    insert_partition_pp(partitions, pp_count, &head_pp, start_index + pp_count, j);
+		    pp_count++;
+		} else if (diff > -(rti_tol_dist)) {
+		    /*
+		     * Make a subtle but important distinction here.  Even
+		     * though the two distances are "equal" within
+		     * tolerance, they are not exactly the same.  If the
+		     * new segment is slightly closer to the ray origin,
+		     * then use its IN point.
+		     *
+		     * This is an attempt to reduce the deflected normals
+		     * sometimes seen along the edges of e.g. a cylinder
+		     * unioned with an ARB8, where the ray hits the top of
+		     * the cylinder and the *side* face of the ARB8 rather
+		     * than the top face of the ARB8.
+		     */
+		    diff = segp->seg_in.hit_dist - pp->inhit.hit_dist;
+		    if (pp->back_pp == head_pp || partitions[pp->back_pp].outhit.hit_dist <=
+			    segp->seg_in.hit_dist) {
+			if (NEAR_ZERO(diff, rti_tol_dist) &&
+				diff < 0) {
+			    pp->inseg = k;
+			    pp->inhit = segp->seg_in;
+			    pp->inflip = 0;
+			}
+		    }
+		} else {
+		    /*
+		     * diff < ~~0
+		     *
+		     * Seg starts before current partition starts,
+		     * but after the previous partition ends.
+		     *	SSSSSSSS...
+		     *	     PPPPP...
+		     *	newpp|pp
+		     */
+		    newpp = &partitions[start_index + pp_count];
+                    initialize_partition(partitions, start_index + pp_count);
+
+		    set(segs_bv, (start_index + pp_count) * bv_index, k-h[id]);
+		    newpp->inseg = lastseg;
+		    newpp->inhit = *lasthit;
+		    newpp->inflip = lastflip;
+		    diff = segp->seg_out.hit_dist - pp->inhit.hit_dist;
+		    if (diff < -(rti_tol_dist)) {
+			/*
+			 * diff < ~0
+			 * Seg starts and ends before current
+			 * partition, but after previous
+			 * partition ends (diff < 0).
+			 *		SSSS
+			 *	pppp		PPPPP...
+			 *		newpp	pp
+			 */
+			newpp->outseg = k;
+			newpp->outhit = segp->seg_out;
+			newpp->outflip = 0;
+			insert_partition_pp(partitions, pp_count, &head_pp, start_index + pp_count, j);
+			pp_count++;
+			break;
+		    } else if (diff < rti_tol_dist) {
+			/*
+			 * diff ~= 0
+			 *
+			 * Seg starts before current
+			 * partition starts, and ends at or
+			 * near the start of the partition.
+			 * (diff == 0).  FUSE the points.
+			 *	SSSSSS
+			 *	     PPPPP
+			 *	newpp|pp
+			 * NOTE: only copy hit point, not
+			 * normals or norm private stuff.
+			 */
+			newpp->outseg = k;
+			newpp->outhit = segp->seg_out;
+			newpp->outhit.hit_dist = pp->inhit.hit_dist;
+			newpp->outflip = 0;
+			insert_partition_pp(partitions, pp_count, &head_pp, start_index + pp_count, j);
+			pp_count++;
+			break;
+		    }
+		    /*
+		     * Seg starts before current partition
+		     * starts, and ends after the start of the
+		     * partition.  (diff > 0).
+		     *	SSSSSSSSSS
+		     *	      PPPPPPP
+		     *	newpp| pp | ...
+		     */
+		    newpp->outseg = pp->inseg;
+		    newpp->outhit = pp->inhit;
+		    newpp->outflip = 1;
+		    lastseg = pp->inseg;
+		    lasthit = &pp->inhit;
+		    lastflip = newpp->outflip;
+		    insert_partition_pp(partitions, pp_count, &head_pp, start_index + pp_count, j);
+		    pp_count++;
+		}
+
+		/*
+		 * Segment and partition start at (roughly) the same
+		 * point.  When fusing 2 points together i.e., when
+		 * NEAR_ZERO(diff, tol) is true, the two points MUST
+		 * be forced to become exactly equal!
+		 */
+		diff = segp->seg_out.hit_dist - pp->outhit.hit_dist;
+		if (diff > rti_tol_dist) {
+		    /*
+		     * Seg & partition start at roughly
+		     * the same spot,
+		     * seg extends beyond partition end.
+		     *	PPPP
+		     *	SSSSSSSS
+		     *	pp  |  newpp
+		     */
+		    set(segs_bv, j * bv_index, k-h[id]);
+
+		    lasthit = &pp->outhit;
+		    lastseg = pp->outseg;
+		    lastflip = 1;
+		    continue;
+		}
+		if (diff > -(rti_tol_dist)) {
+		    /*
+		     * diff ~= 0
+		     * Segment and partition start & end
+		     * (nearly) together.
+		     *	 PPPP
+		     *	 SSSS
+		     */
+		    set(segs_bv, j * bv_index, k-h[id]);
+		    break;
+		} else {
+		    /*
+		     * diff < ~0
+		     *
+		     * Segment + Partition start together,
+		     * segment ends before partition ends.
+		     *	PPPPPPPPPP
+		     *	SSSSSS
+		     *	newpp| pp
+		     */
+		    newpp = &partitions[start_index + pp_count];
+		    *newpp = *pp;
+		    copy_bv(segs_bv, bv_index, (start_index + pp_count) * bv_index, j * bv_index);
+
+		    set(segs_bv, (start_index + pp_count) * bv_index, k-h[id]);
+		    newpp->outseg = k;
+		    newpp->outhit = segp->seg_out;
+		    newpp->outflip = 0;
+		    pp->inseg = k;
+		    pp->inhit = segp->seg_out;
+		    pp->inflip = 1;
+		    insert_partition_pp(partitions, pp_count, &head_pp, start_index + pp_count, j);
+		    pp_count++;
+		    break;
+		}
+		/* NOTREACHED */
+	    }
+	    /*
+	     * Segment has portion which extends beyond the end
+	     * of the last partition.  Tack on the remainder.
+	     *  	PPPPP
+	     *  	     SSSSS
+	     */
+	    if (pp_count > 0 && j == UINT_MAX) {
+		newpp = &partitions[start_index + pp_count];
+                initialize_partition(partitions, start_index + pp_count);
+
+		set(segs_bv, (start_index + pp_count) * bv_index, k-h[id]);
+		newpp->inseg = lastseg;
+		newpp->inhit = *lasthit;
+		newpp->inflip = lastflip;
+		newpp->outseg = k;
+		newpp->outhit = segp->seg_out;
+		append_partition_pp(partitions, pp_count, start_index + pp_count, &tail_pp);
+		pp_count++;
+	    }
+	}
+    }
+
+    if (pp_count > 0) {
+	/* Store the head index of the first partition in this ray */
+	head_partition[id] = head_pp;
+    }
+}
+
+int
+bool_eval(global struct partition *partitions, RESULT_TYPE segs, global uint *h,
+        global uint *segs_bv, const uint bv_index, uint offset, size_t id,
+        global struct bool_region *bregions, global struct tree_bit *btree, const uint region_index)
+{
+    int sp[BOOL_STACKSIZE];
+    int ret;
+    int stackend;
+    uint uop;
+    int idx;
+
+    stackend = 0;
+    sp[stackend++] = INT_MAX;
+    idx = bregions[region_index].btree_offset;
+    for(;;) {
+        for (;;) {
+            uop = btree[idx].val & 7;
+
+            switch (uop) {
+                case UOP_SOLID:
+                    {
+                        /* Tree Leaf */
+                        const uint st_bit = btree[idx].val >> 3;
+                        global struct partition *pp;
+                        RESULT_TYPE segp;
+                        ret = 0;
+
+                        pp = &partitions[offset];
+                        /* Iterate over segments of partition */
+                        for (uint i = 0; i < bv_index; i++) {
+                            uint mask = segs_bv[offset * bv_index + i];
+                            while (mask != 0) {
+                                uint lz = clz(mask);
+                                uint k = h[id] + (31 - lz);
+                                if (isset(segs_bv, offset * bv_index, k - h[id]) != 0) {
+                                    segp = segs+k;
+
+                                    if (segp->seg_sti == st_bit) {
+                                        ret = 1;
+                                        break;
+                                    }
+                                }
+                                // clear bit in mask
+                                mask &= ~(1 << (31-lz));
+                            }
+                            if (ret) break;
+                        }
+                    }
+                    break;
+
+                case UOP_UNION:
+                case UOP_INTERSECT:
+                case UOP_SUBTRACT:
+                case UOP_XOR:
+                    sp[stackend++] = idx;
+                    idx++;
+                    continue;
+                default:
+                    /* bad sp op */
+                    return BOOL_TRUE;       /* Screw up output */
+            }
+            break;
+        }
+
+        for (;;) {
+            idx = sp[--stackend];
+
+            switch (idx) {
+                case INT_MAX:
+                    return ret;		/* top of tree again */
+                case -1:
+                    /* Special operation for subtraction */
+		    ret = !ret;
+		    continue;
+                case -2:
+                    /*
+		     * Special operation for XOR.  lhs was true.  If rhs
+		     * subtree was true, an overlap condition exists (both
+		     * sides of the XOR are BOOL_TRUE).  Return error
+		     * condition.  If subtree is false, then return BOOL_TRUE
+		     * (from lhs).
+		     */
+		    if (ret) {
+			/* stacked temp val: rhs */
+			return -1;	/* GUARD error */
+		    }
+		    ret = BOOL_TRUE;
+		    stackend--;			/* pop temp val */
+		    continue;
+                case -3:
+                    /*
+		     * Special NOP for XOR.  lhs was false.  If rhs is true,
+		     * take note of its regionp.
+		     */
+		    stackend--;			/* pop temp val */
+		    continue;
+                default:
+                    break;
+            }
+
+            uop = btree[idx].val & 7;
+
+	    /*
+	     * Here, each operation will look at the operation just completed
+	     * (the left branch of the tree generally), and rewrite the top of
+	     * the stack and/or branch accordingly.
+	     */
+            switch (uop) {
+                case UOP_SOLID:
+                    /* bool_eval:  pop SOLID? */
+                    return BOOL_TRUE;	    /* screw up output */
+                case UOP_UNION:
+                    if (ret) continue;	    /* BOOL_TRUE, we are done */
+                    /* lhs was false, rewrite as rhs tree */
+                    idx = btree[idx].val >> 3;
+                    break;
+		case UOP_INTERSECT:
+		    if (!ret) {
+		        ret = BOOL_FALSE;
+			continue;
+		    }
+		    /* lhs was true, rewrite as rhs tree */
+                    idx = btree[idx].val >> 3;
+		    break;
+                case UOP_SUBTRACT:
+		    if (!ret) continue;	/* BOOL_FALSE, we are done */
+		    /* lhs was true, rewrite as NOT of rhs tree */
+		    /* We introduce the special NOT operator here */
+                    sp[stackend++] = -1;
+                    idx = btree[idx].val >> 3;
+		    break;
+		case UOP_XOR:
+		    if (ret) {
+		        /* lhs was true, rhs better not be, or we have an
+			 * overlap condition.  Rewrite as guard node followed
+			 * by rhs.
+			 */
+                        idx = btree[idx].val >> 3;
+			sp[stackend++] = idx;		/* temp val for guard node */
+			sp[stackend++] = -2;
+		    } else {
+			/* lhs was false, rewrite as xnop node and result of
+			 * rhs.
+			 */
+                        idx = btree[idx].val >> 3;
+			sp[stackend++] = idx;		/* temp val for xnop */
+			sp[stackend++] = -3;
+		    }
+		    break;
+		default:
+		    /* bool_eval:  bad pop op */
+		    return BOOL_TRUE;	    /* screw up output */
+            }
+	    break;
+	}
+    }
+    /* NOTREACHED */
+}
+
+/**
+ * For each segment's solid that lies in this partition, add
+ * the list of regions that refer to that solid into the
+ * "regiontable" bitarray.
+ */
+void
+build_regiontable(global uint *regions_table, RESULT_TYPE segs,
+        global uint *segs_bv, global uint *regiontable, const uint pp_idx, const uint seg_idx,
+        const uint bv_index, const uint rt_index, const size_t id)
+{
+    RESULT_TYPE segp;
+
+    /* Iterate over segments of partition */
+    for (uint i = 0; i < bv_index; i++) {
+	uint mask = segs_bv[pp_idx * bv_index + i];
+	while (mask != 0) {
+            uint lz = clz(mask);
+	    uint k = seg_idx + (31 - lz);
+	    if (isset(segs_bv, pp_idx * bv_index, k - seg_idx) != 0) {
+		segp = segs+k;
+
+		/* Search for all regions involved in this partition */
+		for (uint m = 0; m < rt_index; m++) {
+		    regiontable[id * rt_index + m] |= regions_table[segp->seg_sti * rt_index + m];
+		}
+	    }
+	    // clear bit in mask
+	    mask &= ~(1 << (31 - lz));
+	}
+    }
+}
+
+void
+reset_regiontable(global uint *regiontable, const size_t id, const uint rt_index)
+{
+    for (uint k = 0; k < rt_index; k++) {
+        regiontable[id * rt_index + k] = 0;
+    }
+}
+
+int
+rt_defoverlap(global struct partition *partitions, const uint pp_idx, global struct bool_region *reg1,
+	      const uint reg1_id, global struct bool_region *reg2, const uint reg2_id, const uint headpp_idx)
+{
+
+    global struct partition *pp;
+
+    /*
+     * Apply heuristics as to which region should claim partition.
+     */
+    if (reg1->reg_aircode != 0) {
+	/* reg1 was air, replace with reg2 */
+	return 2;
+    }
+
+    pp = &partitions[pp_idx];
+    if (pp->back_pp != headpp_idx) {
+	if (partitions[pp->back_pp].region_id == reg1_id)
+	    return 1;
+	if (partitions[pp->back_pp].region_id == reg2_id)
+	    return 2;
+    }
+
+    /* To provide some consistency from ray to ray, use lowest bit # */
+    if (reg1->reg_bit < reg2->reg_bit)
+	return 1;
+    return 2;
+}
+
+void
+rt_default_multioverlap(global struct partition *partitions, global struct bool_region *bregions, global uint *regiontable,
+			const uint first_region, const uint pp_idx, const uint total_regions, const uint headpp_idx, const size_t id)
+{
+    global struct bool_region *lastregion;
+    int code;
+
+    uint rt_index = total_regions/32 +1;
+    uint lastregion_id;
+
+    // Get first region of the regiontable
+    lastregion = &bregions[first_region];
+    lastregion_id = first_region;
+
+    /* Examine the overlapping regions, pairwise */
+    for (uint i = 0; i < rt_index; i++) {
+	uint mask = regiontable[id * rt_index + i];
+	while (mask != 0) {
+	    uint lz = clz(mask);
+	    uint k = (i * 32) + (31 - lz);
+	    if (k != lastregion_id && isset(regiontable, id * rt_index, k) != 0) {
+		global struct bool_region *regp;
+		regp = &bregions[k];
+
+		code = -1;
+		/*
+		 * Two or more regions claim this partition
+		 */
+		if (lastregion->reg_aircode != 0 && regp->reg_aircode == 0) {
+		    /* last region is air, replace with solid regp */
+		    code = 2;
+		} else if (lastregion->reg_aircode == 0 && regp->reg_aircode != 0) {
+		    /* last region solid, regp is air, keep last */
+		    code = 1;
+		} else if (lastregion->reg_aircode != 0 &&
+			   regp->reg_aircode != 0 &&
+			   regp->reg_aircode == lastregion->reg_aircode) {
+		    /* both are same air, keep last */
+		    code = 1;
+		} else {
+		    /*
+		     * Hand overlap to old-style application-specific
+		     * overlap handler, or default.
+		     * 0 = destroy partition,
+		     * 1 = keep part, claiming region=lastregion
+		     * 2 = keep part, claiming region=regp
+		     */
+		    code = rt_defoverlap(partitions, pp_idx, lastregion, lastregion_id, regp, k, headpp_idx);
+		}
+
+		/* Implement the policy in "code" */
+		switch (code) {
+		    case 0:
+			/*
+			 * Destroy the whole partition.
+			 * Reset regiontable
+			 */
+			reset_regiontable(regiontable, id, rt_index);
+			return;
+		    case 1:
+			/* Keep partition, claiming region = lastregion */
+			clr(regiontable, id * rt_index, k);
+			break;
+		    case 2:
+			/* Keep partition, claiming region = regp */
+			clr(regiontable, id * rt_index, lastregion_id);
+			lastregion = regp;
+			lastregion_id = k;
+			break;
+		}
+	    }
+	    // clear bit in mask
+	    mask &= ~(1 << (31 - lz));
+	}
+    }
+}
+
+__kernel void
+rt_boolfinal(global struct partition *partitions, global uint *head_partition, RESULT_TYPE segs,
+        global uint *h, global uint *segs_bv, const int max_depth,
+        global struct bool_region *bregions, const uint total_regions, global struct tree_bit *rtree,
+        global uint *regiontable, const int cur_pixel, const int last_pixel,
+        global uint *regions_table)
+{
+    const size_t id = get_global_size(0)*get_global_id(1)+get_global_id(0);
+
+    if (id >= (last_pixel-cur_pixel+1))
+	return;
+
+    uint head;
+    uint lastregion_idx;
+    int claiming_regions;
+    double diff;
+    uint bv_index = max_depth/32 + 1;
+    uint rt_index = total_regions/32 +1;
+
+    uint lastpp_eval_idx = UINT_MAX;
+
+    //No partitions
+    if (head_partition[id] == UINT_MAX) {
+	return;
+    }
+
+    //Get first partition of the ray
+    head = head_partition[id];
+    head_partition[id] = UINT_MAX;
+
+    //iterate over partitions
+    for (uint current_index = head; current_index != UINT_MAX; current_index = partitions[current_index].forw_pp) {
+	global struct partition *pp = &partitions[current_index];
+        uint first_region_idx = UINT_MAX;
+
+	claiming_regions = 0;
+	/* Force "very thin" partitions to have exactly zero thickness. */
+	if (NEAR_EQUAL(pp->inhit.hit_dist, pp->outhit.hit_dist, rti_tol_dist)) {
+	    pp->outhit.hit_dist = pp->inhit.hit_dist;
+	}
+
+	/* Sanity checks on sorting. */
+	if (pp->inhit.hit_dist > pp->outhit.hit_dist) {
+	    /* Inverted partition */
+	    return;
+	}
+
+	if (pp->forw_pp != UINT_MAX) {
+	    diff = pp->outhit.hit_dist - partitions[pp->forw_pp].inhit.hit_dist;
+	    if (!ZERO(diff)) {
+		if (NEAR_ZERO(diff, rti_tol_dist)) {
+		    /* Fusing 2 partitions */
+		    partitions[pp->forw_pp].inhit.hit_dist = pp->outhit.hit_dist;
+		} else if (diff > 0) {
+		    /* Sorting defect */
+		    return;
+		}
+	    }
+	}
+
+	/* Start with a clean state when evaluating this partition */
+	reset_regiontable(regiontable, id, rt_index);
+
+	/*
+	 * Build regiontable bitarray of all the regions involved in this
+	 * partitions to later evaluate the partitions against the involved
+	 * regions and to resolve any overlap that may occur
+	 */
+	build_regiontable(regions_table, segs, segs_bv, regiontable, current_index, h[id], bv_index, rt_index, id);
+
+	/* Evaluate the boolean trees of any regions involved */
+	for (uint i = 0; i < rt_index; i++) {
+	    uint mask = regiontable[id * rt_index + i];
+	    while (mask != 0) {
+		uint lz = clz(mask);
+		uint k = (i * 32) + (31 - lz);
+		if (isset(regiontable, id * rt_index, k) != 0) {
+
+                    if (first_region_idx == UINT_MAX)
+                        first_region_idx = k;
+
+		    if (bregions[k].reg_all_unions) {
+			claiming_regions++;
+			lastregion_idx = k;
+		    }
+
+		    else if (bool_eval(partitions, segs, h, segs_bv, bv_index, current_index, id, bregions, rtree, k) == BOOL_TRUE) {
+			/* This region claims partition */
+			claiming_regions++;
+			lastregion_idx = k;
+		    }
+		}
+		// clear bit in mask
+		mask &= ~(1 << (31 - lz));
+	    }
+	}
+
+	if (claiming_regions == 0)
+	    continue;
+
+	if (claiming_regions > 1) {
+	    /* There is an overlap between two or more regions */
+	    rt_default_multioverlap(partitions, bregions, regiontable, first_region_idx, current_index, total_regions, head, id);
+
+	    /* Count number of remaining regions, s/b 0 or 1 */
+	    claiming_regions = 0;
+	    for (uint i = 0; i < rt_index; i++) {
+		uint mask = regiontable[id * rt_index + i];
+		while (mask != 0) {
+		    uint lz = clz(mask);
+		    uint k = (i * 32) + (31 - lz);
+
+		    if (isset(regiontable, id * rt_index, k) != 0) {
+			claiming_regions++;
+			lastregion_idx = k;
+		    }
+		    // clear bit in mask
+		    mask &= ~(1 << (31 - lz));
+		}
+	    }
+
+	    /* If claiming_regions != 1, discard partition. */
+	    if (claiming_regions != 1)
+		continue;
+	}
+
+	/*
+	 * claiming_regions == 1
+	 *
+	 * Partition evaluated
+	 */
+	{
+	    global struct partition *lastpp;
+
+            if (head_partition[id] == UINT_MAX) {
+                /* First partition evaluated for this ray
+                 * Start shading at this partition index
+                 */
+                head_partition[id] = current_index;
+            }
+
+	    /* Record the "owning" region. */
+	    pp->region_id = lastregion_idx;
+
+	    /* See if this new partition extends the previous last
+	     * partition, "exactly" matching.
+	     */
+	    if (lastpp_eval_idx != UINT_MAX) {
+		/* there is one last partition evaluated for this ray */
+		lastpp = &partitions[lastpp_eval_idx];
+                lastpp->forw_pp = current_index;
+	    }
+
+	    if (lastpp_eval_idx != UINT_MAX && lastregion_idx == lastpp->region_id &&
+		    NEAR_EQUAL(pp->inhit.hit_dist, lastpp->outhit.hit_dist, rti_tol_dist)) {
+		/* same region, merge by extending last final partition */
+		lastpp->outhit = pp->outhit;
+		lastpp->outseg = pp->outseg;
+		lastpp->outflip = pp->outflip;
+
+		/* Don't lose the fact that the two solids of this
+		 * partition contributed.
+		 */
+		set(segs_bv, lastpp_eval_idx + (bv_index - 1), pp->inseg - h[id]);
+		set(segs_bv, lastpp_eval_idx + (bv_index - 1), pp->outseg - h[id]);
+	    } else {
+		lastpp_eval_idx = current_index;
+	    }
+	}
+        return;
+    }
+}
+
+#endif
+
+/*
+ * Local Variables:
+ * mode: C
+ * tab-width: 8
+ * indent-tabs-mode: t
+ * c-file-style: "stroustrup"
+ * End:
+ * ex: shiftwidth=4 tabstop=8
+ */
+
Index: src/librt/primitives/common.cl
===================================================================
--- src/librt/primitives/common.cl	(revision 70095)
+++ src/librt/primitives/common.cl	(working copy)
@@ -57,6 +57,40 @@
     uint seg_sti;
 };
 
+struct partition {
+    struct hit inhit;
+    struct hit outhit;
+    uint inseg;
+    uint outseg;
+    uint forw_pp;               /* index to the next partition */
+    uint back_pp;               /* index to the previous partition */
+    uint region_id;             /* id of the "owning" region */
+    char inflip;		/* flip inhit->hit_normal */
+    char outflip;		/* flip outhit->hit_normal */
+};
+
+/**
+ * bit expr tree representation
+ *
+ * node:
+ *      uint uop : 3
+ *      uint right_child : 29
+ *
+ * leaf:
+ *      uint uop : 3
+ *      uint st_bit : 29
+ */
+struct tree_bit {
+    uint val;
+};
+
+struct bool_region {
+    uint btree_offset;          /* index to the start of the bit tree */
+    int reg_aircode;            /* Region ID AIR code */
+    int reg_bit;                /* constant index into Regions[] */
+    short reg_all_unions;       /* 1=boolean tree is all unions */
+};
+
 struct region;
 #if RT_SINGLE_HIT
 struct accum {
@@ -70,6 +104,7 @@
     global uint *indexes;
     global uchar *prims;
     global struct region *regions;
+    global int *iregions;
     int lightmodel;
 };
 
Index: src/librt/primitives/primitive_util.c
===================================================================
--- src/librt/primitives/primitive_util.c	(revision 70095)
+++ src/librt/primitives/primitive_util.c	(working copy)
@@ -26,9 +26,7 @@
 
 #include "bu/malloc.h"
 #include "bu/opt.h"
-#include "optical.h"
 #include "../librt_private.h"
-#include "optical/plastic.h"
 
 
 /**
@@ -477,6 +475,7 @@
 static cl_program clt_sh_program, clt_mh_program;
 static cl_kernel clt_frame_kernel;
 static cl_kernel clt_count_hits_kernel, clt_store_segs_kernel, clt_shade_segs_kernel;
+static cl_kernel clt_boolweave_kernel, clt_boolfinal_kernel;
 
 static size_t max_wg_size;
 static cl_uint max_compute_units;
@@ -483,8 +482,10 @@
 
 static cl_mem clt_rand_halftab;
 
-static cl_mem clt_db_ids, clt_db_indexes, clt_db_prims, clt_db_bvh, clt_db_regions;
+static cl_mem clt_db_ids, clt_db_indexes, clt_db_prims, clt_db_bvh, clt_db_regions, clt_db_iregions;
+static cl_mem clt_db_rtree, clt_db_bool_regions, clt_db_regions_table;
 static cl_uint clt_db_nprims;
+static cl_uint clt_db_nregions;
 
 
 
@@ -594,6 +595,8 @@
 
     clReleaseKernel(clt_count_hits_kernel);
     clReleaseKernel(clt_store_segs_kernel);
+    clReleaseKernel(clt_boolweave_kernel);
+    clReleaseKernel(clt_boolfinal_kernel);
     clReleaseKernel(clt_shade_segs_kernel);
 
     clReleaseKernel(clt_frame_kernel);
@@ -615,6 +618,7 @@
     if (!clt_initialized) {
         const char *main_files[] = {
             "solver.cl",
+            "bool.cl",
 
             "arb8_shot.cl",
             "bot_shot.cl",
@@ -662,6 +666,10 @@
         if (error != CL_SUCCESS) bu_bomb("failed to create an OpenCL kernel");
         clt_store_segs_kernel = clCreateKernel(clt_mh_program, "store_segs", &error);
         if (error != CL_SUCCESS) bu_bomb("failed to create an OpenCL kernel");
+	clt_boolweave_kernel = clCreateKernel(clt_mh_program, "rt_boolweave", &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create an OpenCL kernel");
+	clt_boolfinal_kernel = clCreateKernel(clt_mh_program, "rt_boolfinal", &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create an OpenCL kernel");
         clt_shade_segs_kernel = clCreateKernel(clt_mh_program, "shade_segs", &error);
         if (error != CL_SUCCESS) bu_bomb("failed to create an OpenCL kernel");
 
@@ -699,30 +707,6 @@
 }
 
 
-/*
- * Values for Shader Function ID.
- */
-#define SH_NONE 	0
-#define SH_PHONG	1
-
-
-struct clt_phong_specific {
-    cl_double wgt_specular;
-    cl_double wgt_diffuse;
-    cl_int shine;
-};
-
-/**
- * The region structure.
- */
-struct clt_region {
-    cl_float color[3];		/**< @brief explicit color:  0..1  */
-    cl_int mf_id;
-    union {
-	struct clt_phong_specific phg_spec;
-    }udata;
-};
-
 void
 clt_db_store(size_t count, struct soltab *solids[])
 {
@@ -730,55 +714,20 @@
 
     if (count != 0) {
         cl_uchar *ids;
-        struct clt_region *regions;
+	cl_int *iregions;
         cl_uint *indexes;
         struct bu_pool *pool;
         size_t i;
-        
+
 	ids = (cl_uchar*)bu_calloc(count, sizeof(*ids), "ids");
-	regions = (struct clt_region*)bu_calloc(count, sizeof(*regions), "regions");
+	iregions = (cl_int*)bu_calloc(count, sizeof(*iregions), "iregions");
 	for (i=0; i < count; i++) {
 	    const struct soltab *stp = solids[i];
-	    const struct region *regp;
-	    const struct mfuncs *mfp;
-	    const cl_float unset[3] = {1.0f, 1.0f, 1.0f};
+	    const struct region *regp = (struct region *)BU_PTBL_GET(&stp->st_regions,0);
 
             ids[i] = stp->st_id;
-
-	    VMOVE(regions[i].color, unset);
-	    regions[i].mf_id = SH_PHONG;
-
-	    if (BU_PTBL_LEN(&stp->st_regions) > 0) {
-		regp = (const struct region*)BU_PTBL_GET(&stp->st_regions, BU_PTBL_LEN(&stp->st_regions)-1);
-		RT_CK_REGION(regp);
-
-		if (regp->reg_mater.ma_color_valid) {
-		    VMOVE(regions[i].color, regp->reg_mater.ma_color);
-		    regions[i].mf_id = SH_PHONG;
-		}
-
-		mfp = (const struct mfuncs*)regp->reg_mfuncs;
-		if (mfp) {
-		    if (bu_strcmp(mfp->mf_name, "default") ||
-			    bu_strcmp(mfp->mf_name, "phong") ||
-			    bu_strcmp(mfp->mf_name, "plastic") ||
-			    bu_strcmp(mfp->mf_name, "mirror") ||
-			    bu_strcmp(mfp->mf_name, "glass")) {
-			struct phong_specific *src =
-			    (struct phong_specific*)regp->reg_udata;
-			struct clt_phong_specific *dst =
-			    &regions[i].udata.phg_spec;
-
-			dst->shine = src->shine;
-			dst->wgt_diffuse = src->wgt_diffuse;
-			dst->wgt_specular = src->wgt_specular;
-
-			regions[i].mf_id = SH_PHONG;
-		    } else {
-			bu_log("Unknown OCL shader: %s\n", mfp->mf_name);
-		    }
-		}
-	    }
+	    RT_CK_REGION(regp);
+	    iregions[i] = regp->reg_bit;
 	}
 
 	indexes = (cl_uint*)bu_calloc(count+1, sizeof(*indexes), "indexes");
@@ -792,8 +741,8 @@
             /*bu_log("\t(%ld bytes)\n",size);*/
 	    indexes[i] = indexes[i-1] + size;
 	}
-        bu_log("OCLDB:\t%ld primitives\n\t%.2f KB indexes, %.2f KB ids, %.2f KB prims, %.2f KB regions\n", count,
-		(sizeof(*indexes)*(count+1))/1024.0, (sizeof(*ids)*count)/1024.0, indexes[count]/1024.0, (sizeof(*regions)*count)/1024.0);
+        bu_log("OCLDB:\t%ld primitives\n\t%.2f KB indexes, %.2f KB ids, %.2f KB prims\n", count,
+		(sizeof(*indexes)*(count+1))/1024.0, (sizeof(*ids)*count)/1024.0, indexes[count]/1024.0);
 
 	if (indexes[count] != 0) {
 	    clt_db_prims = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, indexes[count], pool->block, &error);
@@ -801,16 +750,16 @@
 	}
         bu_pool_delete(pool);
 
-	clt_db_ids = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(cl_uchar)*count, ids, &error);
+	clt_db_ids = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(*ids)*count, ids, &error);
 	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL ids buffer");
-	clt_db_regions = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(struct clt_region)*count, regions, &error);
-	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL regions buffer");
-	clt_db_indexes = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(cl_uint)*(count+1), indexes, &error);
+	clt_db_iregions = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(*iregions)*count, iregions, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL iregions buffer");
+	clt_db_indexes = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(*indexes)*(count+1), indexes, &error);
 	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL indexes buffer");
 
 	bu_free(indexes, "indexes");
+	bu_free(iregions, "iregions");
 	bu_free(ids, "ids");
-	bu_free(regions, "regions");
     }
 
     clt_db_nprims = count;
@@ -827,15 +776,51 @@
 }
 
 void
+clt_db_store_regions(size_t sz_btree_array, struct bit_tree *btp, size_t nregions, struct cl_bool_region *regions, struct cl_region *mtls)
+{
+    cl_int error;
+
+    if (nregions != 0) {
+	bu_log("OCLRegions:\t%ld regions\n\t%.2f KB regions, %.2f KB mtls\n", nregions, (sizeof(*regions)*nregions)/1024.0, (sizeof(*mtls)*nregions)/1024.0);
+	clt_db_bool_regions = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(*regions)*nregions, regions, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL boolean regions buffer");
+
+	clt_db_rtree = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(*btp)*sz_btree_array, btp, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL boolean trees buffer");
+
+	clt_db_regions = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(*mtls)*nregions, mtls, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL regions buffer");
+    }
+
+    clt_db_nregions = nregions;
+}
+
+void
+clt_db_store_regions_table(cl_uint *regions_table, size_t regions_table_size)
+{
+    cl_int error;
+
+    if (regions_table_size != 0) {
+	clt_db_regions_table = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(cl_uint)*regions_table_size, regions_table, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL regions_table buffer");
+    }
+}
+
+void
 clt_db_release(void)
 {
+    clReleaseMemObject(clt_db_iregions);
     clReleaseMemObject(clt_db_regions);
     clReleaseMemObject(clt_db_bvh);
     clReleaseMemObject(clt_db_prims);
     clReleaseMemObject(clt_db_indexes);
     clReleaseMemObject(clt_db_ids);
+    clReleaseMemObject(clt_db_rtree);
+    clReleaseMemObject(clt_db_bool_regions);
+    clReleaseMemObject(clt_db_regions_table);
 
     clt_db_nprims = 0;
+    clt_db_nregions = 0;
 }
 
 void
@@ -843,7 +828,7 @@
 	  int width, int ibackground[3], int inonbackground[3],
 	  double airdensity, double haze[3], fastf_t gamma,
           mat_t view2model, fastf_t cell_width, fastf_t cell_height,
-          fastf_t aspect, int lightmodel)
+          fastf_t aspect, int lightmodel, int a_no_booleans)
 {
     const size_t npix = last_pixel-cur_pixel+1;
 
@@ -898,156 +883,231 @@
 
     bu_log("%ldx%ld grid, %ldx%ld subgrids\n", wxh[0], wxh[1], swxh[0], swxh[1]);
 
-    switch (lightmodel) {
-	case 5:
-	    {
-	    size_t sz_counts;
-	    cl_int *counts;
-	    cl_mem pcounts;
-	    size_t sz_h;
-	    cl_uint *h;
-	    cl_mem ph;
-	    size_t sz_segs;
-	    cl_mem psegs;
-	    size_t snpix = swxh[0]*swxh[1];
+    if (a_no_booleans) {
+	bu_semaphore_acquire(clt_semaphore);
+	error = clSetKernelArg(clt_frame_kernel, 0, sizeof(cl_mem), &ppixels);
+	error |= clSetKernelArg(clt_frame_kernel, 1, sizeof(cl_uchar3), &p.o);
+	error |= clSetKernelArg(clt_frame_kernel, 2, sizeof(cl_int), &p.cur_pixel);
+	error |= clSetKernelArg(clt_frame_kernel, 3, sizeof(cl_int), &p.last_pixel);
+	error |= clSetKernelArg(clt_frame_kernel, 4, sizeof(cl_int), &p.width);
+	error |= clSetKernelArg(clt_frame_kernel, 5, sizeof(cl_mem), &clt_rand_halftab);
+	error |= clSetKernelArg(clt_frame_kernel, 6, sizeof(cl_uint), &p.randhalftabsize);
+	error |= clSetKernelArg(clt_frame_kernel, 7, sizeof(cl_uchar3), &p.ibackground);
+	error |= clSetKernelArg(clt_frame_kernel, 8, sizeof(cl_uchar3), &p.inonbackground);
+	error |= clSetKernelArg(clt_frame_kernel, 9, sizeof(cl_double), &p.airdensity);
+	error |= clSetKernelArg(clt_frame_kernel, 10, sizeof(cl_double3), &p.haze);
+	error |= clSetKernelArg(clt_frame_kernel, 11, sizeof(cl_double), &p.gamma);
+	error |= clSetKernelArg(clt_frame_kernel, 12, sizeof(cl_double16), &p.view2model);
+	error |= clSetKernelArg(clt_frame_kernel, 13, sizeof(cl_double), &p.cell_width);
+	error |= clSetKernelArg(clt_frame_kernel, 14, sizeof(cl_double), &p.cell_height);
+	error |= clSetKernelArg(clt_frame_kernel, 15, sizeof(cl_double), &p.aspect);
+	error |= clSetKernelArg(clt_frame_kernel, 16, sizeof(cl_int), &lightmodel);
+	error |= clSetKernelArg(clt_frame_kernel, 17, sizeof(cl_uint), &clt_db_nprims);
+	error |= clSetKernelArg(clt_frame_kernel, 18, sizeof(cl_mem), &clt_db_ids);
+	error |= clSetKernelArg(clt_frame_kernel, 19, sizeof(cl_mem), &clt_db_bvh);
+	error |= clSetKernelArg(clt_frame_kernel, 20, sizeof(cl_mem), &clt_db_indexes);
+	error |= clSetKernelArg(clt_frame_kernel, 21, sizeof(cl_mem), &clt_db_prims);
+	error |= clSetKernelArg(clt_frame_kernel, 22, sizeof(cl_mem), &clt_db_regions);
+	error |= clSetKernelArg(clt_frame_kernel, 23, sizeof(cl_mem), &clt_db_iregions);
+	if (error != CL_SUCCESS) bu_bomb("failed to set OpenCL kernel arguments");
+	error = clEnqueueNDRangeKernel(clt_queue, clt_frame_kernel, 2, NULL, wxh,
+		swxh, 0, NULL, NULL);
+	bu_semaphore_release(clt_semaphore);
+    } else {
+	size_t sz_counts;
+	cl_int *counts;
+	cl_mem pcounts;
+	size_t sz_h;
+	cl_uint *h;
+	cl_mem ph;
+	size_t sz_segs;
+	cl_mem psegs;
+	size_t sz_ipartitions;
+	cl_uint *ipart;
+	cl_mem head_partition;
+	size_t sz_partitions;
+	cl_mem ppartitions;
+	cl_int max_depth;
+	size_t sz_bv;
+	cl_uint *bv;
+	cl_mem segs_bv;
+	size_t sz_regiontable;
+	cl_uint *regiontable;
+	cl_mem regiontable_bv;
+	size_t snpix = swxh[0]*swxh[1];
 
-	    sz_counts = sizeof(cl_int)*npix;
-	    pcounts = clCreateBuffer(clt_context, CL_MEM_WRITE_ONLY|CL_MEM_HOST_READ_ONLY, sz_counts, NULL, &error);
-	    if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL counts buffer");
+	sz_counts = sizeof(cl_int)*npix;
+	pcounts = clCreateBuffer(clt_context, CL_MEM_WRITE_ONLY|CL_MEM_HOST_READ_ONLY, sz_counts, NULL, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL counts buffer");
 
+	bu_semaphore_acquire(clt_semaphore);
+	error = clSetKernelArg(clt_count_hits_kernel, 0, sizeof(cl_mem), &pcounts);
+	error |= clSetKernelArg(clt_count_hits_kernel, 1, sizeof(cl_int), &p.cur_pixel);
+	error |= clSetKernelArg(clt_count_hits_kernel, 2, sizeof(cl_int), &p.last_pixel);
+	error |= clSetKernelArg(clt_count_hits_kernel, 3, sizeof(cl_int), &p.width);
+	error |= clSetKernelArg(clt_count_hits_kernel, 4, sizeof(cl_double16), &p.view2model);
+	error |= clSetKernelArg(clt_count_hits_kernel, 5, sizeof(cl_double), &p.cell_width);
+	error |= clSetKernelArg(clt_count_hits_kernel, 6, sizeof(cl_double), &p.cell_height);
+	error |= clSetKernelArg(clt_count_hits_kernel, 7, sizeof(cl_double), &p.aspect);
+	error |= clSetKernelArg(clt_count_hits_kernel, 8, sizeof(cl_uint), &clt_db_nprims);
+	error |= clSetKernelArg(clt_count_hits_kernel, 9, sizeof(cl_mem), &clt_db_ids);
+	error |= clSetKernelArg(clt_count_hits_kernel, 10, sizeof(cl_mem), &clt_db_bvh);
+	error |= clSetKernelArg(clt_count_hits_kernel, 11, sizeof(cl_mem), &clt_db_indexes);
+	error |= clSetKernelArg(clt_count_hits_kernel, 12, sizeof(cl_mem), &clt_db_prims);
+	if (error != CL_SUCCESS) bu_bomb("failed to set OpenCL kernel arguments");
+	error = clEnqueueNDRangeKernel(clt_queue, clt_count_hits_kernel, 2, NULL, wxh,
+		swxh, 0, NULL, NULL);
+	bu_semaphore_release(clt_semaphore);
+
+	/* once we can do the scan on the device we won't need these transfers */
+	counts = (cl_int*)bu_calloc(1, sz_counts, "counts");
+	clEnqueueReadBuffer(clt_queue, pcounts, CL_TRUE, 0, sz_counts, counts, 0, NULL, NULL);
+	clReleaseMemObject(pcounts);
+
+	sz_h = sizeof(cl_uint)*(npix+1);
+	h = (cl_uint*)bu_calloc(1, sz_h, "h");
+	h[0] = 0;
+	max_depth = 0;
+	for (i=1; i<=npix; i++) {
+	    BU_ASSERT((counts[i-1] % 2) == 0);
+	    h[i] = h[i-1] + counts[i-1]/2;	/* number of segs is half the number of hits */
+	    if (counts[i-1]/2 > max_depth)
+		max_depth = counts[i-1]/2;
+	}
+	bu_free(counts, "counts");
+
+	ph = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sz_h, h, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL offs buffer");
+
+	sz_segs = sizeof(struct cl_seg)*h[npix];
+
+	sz_ipartitions = sizeof(cl_uint)*npix; /* store index to first partition of the ray */
+	sz_partitions = sizeof(struct cl_partition)*h[npix]*2; /*create partition buffer with size= 2*number of segments */
+	ipart = (cl_uint*)bu_calloc(1, sz_ipartitions, "ipart");
+
+	sz_bv = sizeof(cl_uint)*(h[npix]*2)*(max_depth/32 + 1); /* bitarray to represent the segs in each partition */
+	bv = (cl_uint*)bu_calloc(1, sz_bv, "bv");
+
+	sz_regiontable = sizeof(cl_uint)*npix*(clt_db_nregions/32 +1); /* bitarray to represent the regions involved in each partition */
+	regiontable = (cl_uint*)bu_calloc(1, sz_regiontable, "regiontable");
+
+	bu_free(h, "h");
+
+	if (sz_segs != 0) {
+	    psegs = clCreateBuffer(clt_context, CL_MEM_READ_WRITE|CL_MEM_HOST_NO_ACCESS, sz_segs, NULL, &error);
+	    if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL segs buffer");
+
 	    bu_semaphore_acquire(clt_semaphore);
-	    error = clSetKernelArg(clt_count_hits_kernel, 0, sizeof(cl_mem), &pcounts);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 1, sizeof(cl_int), &p.cur_pixel);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 2, sizeof(cl_int), &p.last_pixel);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 3, sizeof(cl_int), &p.width);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 4, sizeof(cl_double16), &p.view2model);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 5, sizeof(cl_double), &p.cell_width);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 6, sizeof(cl_double), &p.cell_height);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 7, sizeof(cl_double), &p.aspect);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 8, sizeof(cl_uint), &clt_db_nprims);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 9, sizeof(cl_mem), &clt_db_ids);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 10, sizeof(cl_mem), &clt_db_bvh);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 11, sizeof(cl_mem), &clt_db_indexes);
-	    error |= clSetKernelArg(clt_count_hits_kernel, 12, sizeof(cl_mem), &clt_db_prims);
+	    error = clSetKernelArg(clt_store_segs_kernel, 0, sizeof(cl_mem), &psegs);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 1, sizeof(cl_mem), &ph);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 2, sizeof(cl_int), &p.cur_pixel);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 3, sizeof(cl_int), &p.last_pixel);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 4, sizeof(cl_int), &p.width);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 5, sizeof(cl_double16), &p.view2model);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 6, sizeof(cl_double), &p.cell_width);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 7, sizeof(cl_double), &p.cell_height);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 8, sizeof(cl_double), &p.aspect);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 9, sizeof(cl_uint), &clt_db_nprims);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 10, sizeof(cl_mem), &clt_db_ids);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 11, sizeof(cl_mem), &clt_db_bvh);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 12, sizeof(cl_mem), &clt_db_indexes);
+	    error |= clSetKernelArg(clt_store_segs_kernel, 13, sizeof(cl_mem), &clt_db_prims);
 	    if (error != CL_SUCCESS) bu_bomb("failed to set OpenCL kernel arguments");
-	    error = clEnqueueNDRangeKernel(clt_queue, clt_count_hits_kernel, 2, NULL, wxh,
+	    error = clEnqueueNDRangeKernel(clt_queue, clt_store_segs_kernel, 2, NULL, wxh,
 		    swxh, 0, NULL, NULL);
 	    bu_semaphore_release(clt_semaphore);
+	} else {
+	    psegs = NULL;
+	}
 
-	    /* once we can do the scan on the device we won't need these transfers */
-	    counts = (cl_int*)bu_calloc(1, sz_counts, "counts");
-	    clEnqueueReadBuffer(clt_queue, pcounts, CL_TRUE, 0, sz_counts, counts, 0, NULL, NULL);
-	    clReleaseMemObject(pcounts);
+	head_partition = clCreateBuffer(clt_context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, sz_ipartitions, ipart, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL head partitions buffer");
+	bu_free(ipart, "ipart");
 
-	    sz_h = sizeof(cl_uint)*(npix+1);
-	    h = (cl_uint*)bu_calloc(1, sz_h, "h");
-	    h[0] = 0;
-	    for (i=1; i<=npix; i++) {
-		BU_ASSERT((counts[i-1] % 2) == 0);
-		h[i] = h[i-1] + counts[i-1]/2;	/* number of segs is half the number of hits */
-	    }
-	    bu_free(counts, "counts");
+	segs_bv = clCreateBuffer(clt_context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, sz_bv, bv, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL segs bitvector buffer");
+	bu_free(bv, "bv");
 
-	    ph = clCreateBuffer(clt_context, CL_MEM_READ_ONLY|CL_MEM_HOST_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sz_h, h, &error);
-	    if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL offs buffer");
+	regiontable_bv = clCreateBuffer(clt_context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR, sz_regiontable, regiontable, &error);
+	if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL segs bitvector buffer");
+	bu_free(regiontable, "regiontable");
 
-	    sz_segs = sizeof(struct cl_seg)*h[npix];
-	    bu_free(h, "h");
+	if (sz_partitions != 0) {
+	    ppartitions = clCreateBuffer(clt_context, CL_MEM_READ_WRITE|CL_MEM_HOST_NO_ACCESS, sz_partitions, NULL, &error);
+	    if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL partitions buffer");
 
-	    if (sz_segs != 0) {
-		psegs = clCreateBuffer(clt_context, CL_MEM_READ_WRITE|CL_MEM_HOST_NO_ACCESS, sz_segs, NULL, &error);
-		if (error != CL_SUCCESS) bu_bomb("failed to create OpenCL segs buffer");
-
-		bu_semaphore_acquire(clt_semaphore);
-		error = clSetKernelArg(clt_store_segs_kernel, 0, sizeof(cl_mem), &psegs);
-		error |= clSetKernelArg(clt_store_segs_kernel, 1, sizeof(cl_mem), &ph);
-		error |= clSetKernelArg(clt_store_segs_kernel, 2, sizeof(cl_int), &p.cur_pixel);
-		error |= clSetKernelArg(clt_store_segs_kernel, 3, sizeof(cl_int), &p.last_pixel);
-		error |= clSetKernelArg(clt_store_segs_kernel, 4, sizeof(cl_int), &p.width);
-		error |= clSetKernelArg(clt_store_segs_kernel, 5, sizeof(cl_double16), &p.view2model);
-		error |= clSetKernelArg(clt_store_segs_kernel, 6, sizeof(cl_double), &p.cell_width);
-		error |= clSetKernelArg(clt_store_segs_kernel, 7, sizeof(cl_double), &p.cell_height);
-		error |= clSetKernelArg(clt_store_segs_kernel, 8, sizeof(cl_double), &p.aspect);
-		error |= clSetKernelArg(clt_store_segs_kernel, 9, sizeof(cl_uint), &clt_db_nprims);
-		error |= clSetKernelArg(clt_store_segs_kernel, 10, sizeof(cl_mem), &clt_db_ids);
-		error |= clSetKernelArg(clt_store_segs_kernel, 11, sizeof(cl_mem), &clt_db_bvh);
-		error |= clSetKernelArg(clt_store_segs_kernel, 12, sizeof(cl_mem), &clt_db_indexes);
-		error |= clSetKernelArg(clt_store_segs_kernel, 13, sizeof(cl_mem), &clt_db_prims);
-		if (error != CL_SUCCESS) bu_bomb("failed to set OpenCL kernel arguments");
-		error = clEnqueueNDRangeKernel(clt_queue, clt_store_segs_kernel, 2, NULL, wxh,
-			swxh, 0, NULL, NULL);
-		bu_semaphore_release(clt_semaphore);
-            } else {
-		psegs = NULL;
-            }
-
 	    bu_semaphore_acquire(clt_semaphore);
-	    error = clSetKernelArg(clt_shade_segs_kernel, 0, sizeof(cl_mem), &ppixels);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 1, sizeof(cl_uchar3), &p.o);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 2, sizeof(cl_mem), &psegs);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 3, sizeof(cl_mem), &ph);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 4, sizeof(cl_int), &p.cur_pixel);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 5, sizeof(cl_int), &p.last_pixel);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 6, sizeof(cl_int), &p.width);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 7, sizeof(cl_mem), &clt_rand_halftab);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 8, sizeof(cl_uint), &p.randhalftabsize);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 9, sizeof(cl_uchar3), &p.ibackground);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 10, sizeof(cl_uchar3), &p.inonbackground);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 11, sizeof(cl_double), &p.airdensity);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 12, sizeof(cl_double3), &p.haze);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 13, sizeof(cl_double), &p.gamma);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 14, sizeof(cl_double16), &p.view2model);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 15, sizeof(cl_double), &p.cell_width);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 16, sizeof(cl_double), &p.cell_height);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 17, sizeof(cl_double), &p.aspect);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 18, sizeof(cl_int), &lightmodel);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 19, sizeof(cl_uint), &clt_db_nprims);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 20, sizeof(cl_mem), &clt_db_ids);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 21, sizeof(cl_mem), &clt_db_bvh);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 22, sizeof(cl_mem), &clt_db_indexes);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 23, sizeof(cl_mem), &clt_db_prims);
-	    error |= clSetKernelArg(clt_shade_segs_kernel, 24, sizeof(cl_mem), &clt_db_regions);
+	    error = clSetKernelArg(clt_boolweave_kernel, 0, sizeof(cl_mem), &ppartitions);
+	    error |= clSetKernelArg(clt_boolweave_kernel, 1, sizeof(cl_mem), &head_partition);
+	    error |= clSetKernelArg(clt_boolweave_kernel, 2, sizeof(cl_mem), &psegs);
+	    error |= clSetKernelArg(clt_boolweave_kernel, 3, sizeof(cl_mem), &ph);
+	    error |= clSetKernelArg(clt_boolweave_kernel, 4, sizeof(cl_mem), &segs_bv);
+	    error |= clSetKernelArg(clt_boolweave_kernel, 5, sizeof(cl_int), &p.cur_pixel);
+	    error |= clSetKernelArg(clt_boolweave_kernel, 6, sizeof(cl_int), &p.last_pixel);
+	    error |= clSetKernelArg(clt_boolweave_kernel, 7, sizeof(cl_int), &max_depth);
 	    if (error != CL_SUCCESS) bu_bomb("failed to set OpenCL kernel arguments");
-	    error = clEnqueueNDRangeKernel(clt_queue, clt_shade_segs_kernel, 1, NULL, &npix,
+	    error = clEnqueueNDRangeKernel(clt_queue, clt_boolweave_kernel, 1, NULL, &npix,
 		    &snpix, 0, NULL, NULL);
 	    bu_semaphore_release(clt_semaphore);
+	} else {
+	    ppartitions = NULL;
+	}
 
-	    clReleaseMemObject(ph);
-	    clReleaseMemObject(psegs);
-	    }
-	    break;
-	default:
-	    {
-	    bu_semaphore_acquire(clt_semaphore);
-	    error = clSetKernelArg(clt_frame_kernel, 0, sizeof(cl_mem), &ppixels);
-	    error |= clSetKernelArg(clt_frame_kernel, 1, sizeof(cl_uchar3), &p.o);
-	    error |= clSetKernelArg(clt_frame_kernel, 2, sizeof(cl_int), &p.cur_pixel);
-	    error |= clSetKernelArg(clt_frame_kernel, 3, sizeof(cl_int), &p.last_pixel);
-	    error |= clSetKernelArg(clt_frame_kernel, 4, sizeof(cl_int), &p.width);
-	    error |= clSetKernelArg(clt_frame_kernel, 5, sizeof(cl_mem), &clt_rand_halftab);
-	    error |= clSetKernelArg(clt_frame_kernel, 6, sizeof(cl_uint), &p.randhalftabsize);
-	    error |= clSetKernelArg(clt_frame_kernel, 7, sizeof(cl_uchar3), &p.ibackground);
-	    error |= clSetKernelArg(clt_frame_kernel, 8, sizeof(cl_uchar3), &p.inonbackground);
-	    error |= clSetKernelArg(clt_frame_kernel, 9, sizeof(cl_double), &p.airdensity);
-	    error |= clSetKernelArg(clt_frame_kernel, 10, sizeof(cl_double3), &p.haze);
-	    error |= clSetKernelArg(clt_frame_kernel, 11, sizeof(cl_double), &p.gamma);
-	    error |= clSetKernelArg(clt_frame_kernel, 12, sizeof(cl_double16), &p.view2model);
-	    error |= clSetKernelArg(clt_frame_kernel, 13, sizeof(cl_double), &p.cell_width);
-	    error |= clSetKernelArg(clt_frame_kernel, 14, sizeof(cl_double), &p.cell_height);
-	    error |= clSetKernelArg(clt_frame_kernel, 15, sizeof(cl_double), &p.aspect);
-	    error |= clSetKernelArg(clt_frame_kernel, 16, sizeof(cl_int), &lightmodel);
-	    error |= clSetKernelArg(clt_frame_kernel, 17, sizeof(cl_uint), &clt_db_nprims);
-	    error |= clSetKernelArg(clt_frame_kernel, 18, sizeof(cl_mem), &clt_db_ids);
-	    error |= clSetKernelArg(clt_frame_kernel, 19, sizeof(cl_mem), &clt_db_bvh);
-	    error |= clSetKernelArg(clt_frame_kernel, 20, sizeof(cl_mem), &clt_db_indexes);
-	    error |= clSetKernelArg(clt_frame_kernel, 21, sizeof(cl_mem), &clt_db_prims);
-	    error |= clSetKernelArg(clt_frame_kernel, 22, sizeof(cl_mem), &clt_db_regions);
-	    if (error != CL_SUCCESS) bu_bomb("failed to set OpenCL kernel arguments");
-	    error = clEnqueueNDRangeKernel(clt_queue, clt_frame_kernel, 2, NULL, wxh,
-		    swxh, 0, NULL, NULL);
-	    bu_semaphore_release(clt_semaphore);
-	    }
-	    break;
+	bu_semaphore_acquire(clt_semaphore);
+	error = clSetKernelArg(clt_boolfinal_kernel, 0, sizeof(cl_mem), &ppartitions);
+        error |= clSetKernelArg(clt_boolfinal_kernel, 1, sizeof(cl_mem), &head_partition);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 2, sizeof(cl_mem), &psegs);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 3, sizeof(cl_mem), &ph);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 4, sizeof(cl_mem), &segs_bv);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 5, sizeof(cl_int), &max_depth);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 6, sizeof(cl_mem), &clt_db_bool_regions);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 7, sizeof(cl_uint), &clt_db_nregions);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 8, sizeof(cl_mem), &clt_db_rtree);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 9, sizeof(cl_mem), &regiontable_bv);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 10, sizeof(cl_int), &p.cur_pixel);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 11, sizeof(cl_int), &p.last_pixel);
+	error |= clSetKernelArg(clt_boolfinal_kernel, 12, sizeof(cl_mem), &clt_db_regions_table);
+	if (error != CL_SUCCESS) bu_bomb("failed to set OpenCL kernel arguments");
+	error = clEnqueueNDRangeKernel(clt_queue, clt_boolfinal_kernel, 1, NULL, &npix,
+		&snpix, 0, NULL, NULL);
+	bu_semaphore_release(clt_semaphore);
+
+	bu_semaphore_acquire(clt_semaphore);
+	error = clSetKernelArg(clt_shade_segs_kernel, 0, sizeof(cl_mem), &ppixels);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 1, sizeof(cl_uchar3), &p.o);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 2, sizeof(cl_mem), &psegs);
+        error |= clSetKernelArg(clt_shade_segs_kernel, 3, sizeof(cl_mem), &head_partition);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 4, sizeof(cl_int), &p.cur_pixel);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 5, sizeof(cl_int), &p.last_pixel);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 6, sizeof(cl_int), &p.width);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 7, sizeof(cl_mem), &clt_rand_halftab);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 8, sizeof(cl_uint), &p.randhalftabsize);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 9, sizeof(cl_uchar3), &p.ibackground);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 10, sizeof(cl_uchar3), &p.inonbackground);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 11, sizeof(cl_double), &p.airdensity);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 12, sizeof(cl_double3), &p.haze);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 13, sizeof(cl_double), &p.gamma);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 14, sizeof(cl_double16), &p.view2model);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 15, sizeof(cl_double), &p.cell_width);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 16, sizeof(cl_double), &p.cell_height);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 17, sizeof(cl_double), &p.aspect);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 18, sizeof(cl_int), &lightmodel);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 19, sizeof(cl_mem), &clt_db_ids);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 20, sizeof(cl_mem), &clt_db_indexes);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 21, sizeof(cl_mem), &clt_db_prims);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 22, sizeof(cl_mem), &clt_db_regions);
+	error |= clSetKernelArg(clt_shade_segs_kernel, 23, sizeof(cl_mem), &ppartitions);
+	if (error != CL_SUCCESS) bu_bomb("failed to set OpenCL kernel arguments");
+	error = clEnqueueNDRangeKernel(clt_queue, clt_shade_segs_kernel, 1, NULL, &npix,
+		&snpix, 0, NULL, NULL);
+	bu_semaphore_release(clt_semaphore);
+
+	clReleaseMemObject(ph);
+	clReleaseMemObject(psegs);
+	clReleaseMemObject(ppartitions);
+	clReleaseMemObject(segs_bv);
+	clReleaseMemObject(regiontable_bv);
+	clReleaseMemObject(head_partition);
     }
     clEnqueueReadBuffer(clt_queue, ppixels, CL_TRUE, 0, sz_pixels, pixels, 0, NULL, NULL);
     clReleaseMemObject(ppixels);
Index: src/librt/primitives/rt.cl
===================================================================
--- src/librt/primitives/rt.cl	(revision 70095)
+++ src/librt/primitives/rt.cl	(working copy)
@@ -317,9 +317,7 @@
 }
 
 inline double3
-shade(const double3 r_pt, const double3 r_dir, struct hit *hitp, const uint idx, const double3 lt_pos,
-      global uchar *ids, global uint *indexes, global uchar *prims, global struct region *regions)
-
+shade(const double3 r_pt, const double3 r_dir, struct hit *hitp, const double3 lt_pos, const uint idx, global struct region *regions)
 {
     double3 color;
     double3 normal;
@@ -364,8 +362,7 @@
     if (acc->lightmodel == 4) {
 	if (seg_in->hit_dist >= 0.0) {
 	    norm(seg_in, acc->r_pt, acc->r_dir, acc->ids[idx], acc->prims + acc->indexes[idx]);
-	    const double3 color = shade(acc->r_pt, acc->r_dir, seg_in, idx,
-		    acc->lt_pos, acc->ids, acc->indexes, acc->prims, acc->regions);
+	    const double3 color = shade(acc->r_pt, acc->r_dir, seg_in, acc->lt_pos, acc->iregions[idx], acc->regions);
 	    double f = exp(-seg_in->hit_dist*1e-10);
 	    acc->a_color += color * f;
 	    acc->a_total += f;
@@ -382,7 +379,7 @@
 	 const double16 view2model, const double cell_width, const double cell_height,
 	 const double aspect, const int lightmodel, const uint nprims, global uchar *ids,
 	 global struct linear_bvh_node *nodes, global uint *indexes,global uchar *prims,
-	 global struct region *regions)
+	 global struct region *regions, global int *iregions)
 {
     const size_t id = get_global_size(0)*get_global_id(1)+get_global_id(0);
 
@@ -412,6 +409,7 @@
     a.indexes = indexes;
     a.prims = prims;
     a.regions = regions;
+    a.iregions = iregions;
 
     struct seg *segp = &a.s;
     segp->seg_in.hit_dist = INFINITY;
@@ -448,7 +446,7 @@
          */
         switch (lightmodel) {
 	    default:
-		a_color = shade(r_pt, r_dir, hitp, idx, lt_pos, ids, indexes, prims, regions);
+		a_color = shade(r_pt, r_dir, hitp, lt_pos, iregions[idx], regions);
 	    	break;
             case 1:
                 /* Light from the "eye" (ray source).  Note sign change */
@@ -596,15 +594,15 @@
 }
 
 __kernel void
-shade_segs(global uchar *pixels, const uchar3 o, RESULT_TYPE segs, global uint *h,
+shade_segs(global uchar *pixels, const uchar3 o, RESULT_TYPE segs, global uint *head_partition,
          const int cur_pixel, const int last_pixel, const int width,
 	 global float *rand_halftab, const uint randhalftabsize,
 	 const uchar3 background, const uchar3 nonbackground,
 	 const double airdensity, const double3 haze, const double gamma,
 	 const double16 view2model, const double cell_width, const double cell_height,
-	 const double aspect, const int lightmodel, const uint nprims, global uchar *ids,
-	 global struct linear_bvh_node *nodes, global uint *indexes, global uchar *prims,
-	 global struct region *regions)
+	 const double aspect, const int lightmodel, global uchar *ids,
+	 global uint *indexes, global uchar *prims, global struct region *regions,
+         global struct partition *partitions)
 {
     const size_t id = get_global_size(0)*get_global_id(1)+get_global_id(0);
 
@@ -626,28 +624,44 @@
     double3 a_color;
     uchar3 rgb;
     struct hit hitp;
+    uint head;
+    bool flipflag;
+    uint region_id;
 
     a_color = 0.0;
     hitp.hit_dist = INFINITY;
-    if (h[id] != h[id+1]) {
+    region_id = 0;
+    if (head_partition[id] != UINT_MAX) {
 	uint idx;
+	double diffuse0 = 0;
+	double cosI0 = 0;
+	double3 work0, work1;
 
 	idx = UINT_MAX;
-	for (uint k=h[id]; k!=h[id+1]; k++) {
-	    RESULT_TYPE segp = segs+k;
 
-	    if (segp->seg_in.hit_dist < hitp.hit_dist) {
-		hitp = segp->seg_in;
-		idx = segp->seg_sti;
-	    }
+	/* Get first partition of the ray */
+	head = head_partition[id];
+	flipflag = 0;
+	for (uint index = head; index != UINT_MAX; index = partitions[index].forw_pp) {
+	    global struct partition *pp = &partitions[index];
+            RESULT_TYPE segp = &segs[pp->inseg];
+
+            if (pp->inhit.hit_dist < hitp.hit_dist) {
+                hitp = pp->inhit;
+                idx = segp->seg_sti;
+                region_id = pp->region_id;
+                flipflag = pp->inflip;
+            }
 	}
+
         double3 normal;
 
-	if (hitp.hit_dist < 0.0) {
+        if (hitp.hit_dist < 0.0) {
 	    /* Eye inside solid, orthoview */
 	    normal = -r_dir;
         } else {
 	    norm(&hitp, r_pt, r_dir, ids[idx], prims + indexes[idx]);
+            hitp.hit_normal = flipflag ? -hitp.hit_normal : hitp.hit_normal;
 	    normal = hitp.hit_normal;
         }
 
@@ -654,8 +668,28 @@
         /*
          * Diffuse reflectance from each light source
          */
-	a_color = shade(r_pt, r_dir, &hitp, idx, lt_pos, ids, indexes, prims, regions);
+        switch (lightmodel) {
+	    default:
+		a_color = shade(r_pt, r_dir, &hitp, lt_pos, region_id, regions);
+	    	break;
+            case 1:
+                /* Light from the "eye" (ray source).  Note sign change */
+                diffuse0 = 0;
+                if ((cosI0 = -dot(normal, r_dir)) >= 0.0)
+                    diffuse0 = cosI0 * (1.0 - AmbientIntensity);
+                work0 = lt_color * diffuse0;
 
+                /* Add in contribution from ambient light */
+                work1 = ambient_color * AmbientIntensity;
+                a_color = work0 + work1;
+                break;
+            case 2:
+                /* Store surface normals pointing inwards */
+                /* (For Spencer's moving light program) */
+                a_color = (normal * DOUBLE_C(-.5)) + DOUBLE_C(.5);
+                break;
+        }
+
         /*
          * e ^(-density * distance)
          */
Index: src/librt/tree.c
===================================================================
--- src/librt/tree.c	(revision 70095)
+++ src/librt/tree.c	(working copy)
@@ -1092,6 +1092,51 @@
 }
 
 
+#ifdef USE_OPENCL
+void
+rt_bit_tree(struct bit_tree *btp, const union tree *tp, size_t *len)
+{
+    int idx;
+    uint st_bit, uop, rchild;
+
+    if (tp == TREE_NULL)
+        return;
+
+    idx = (*len)++;
+    switch (tp->tr_op) {
+        case OP_SOLID:
+            /* Tree Leaf */
+            st_bit = tp->tr_a.tu_stp->st_bit;
+            if (btp) btp[idx].val = (st_bit << 3) | UOP_SOLID;
+            return;
+        case OP_SUBTRACT:
+            uop = UOP_SUBTRACT;
+            break;
+        case OP_UNION:
+            uop = UOP_UNION;
+            break;
+        case OP_INTERSECT:
+            uop = UOP_INTERSECT;
+            break;
+        case OP_XOR:
+            uop = UOP_XOR;
+            break;
+        default:
+            bu_log("rt_bit_tree: bad op[%d]\n", tp->tr_op);
+            exit(1);
+            break;
+    }
+
+    rt_bit_tree(btp, tp->tr_b.tb_left, len);
+
+    rchild = *len;
+    if (btp) btp[idx].val = (rchild << 3) | uop;
+
+    rt_bit_tree(btp, tp->tr_b.tb_right, len);
+}
+#endif
+
+
 /*
  * Local Variables:
  * mode: C

