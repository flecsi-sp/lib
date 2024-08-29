.. |br| raw:: html

   <br />

.. _orientation:

Orientation
***********

Consider the quadrilateral zones :math:`\alpha \colon \{a,b,e,f\}` and
:math:`\beta \colon \{b,c,d,e\}`:

.. tikz::
  :include: orientation/quads.tikz
  :xscale: 50
  :align: center

Depending on the global ordering of :math:`\alpha` and :math:`\beta`,
the handedness of the *interface* (edge) :math:`\overrightarrow{be}`
will result in a positive or negative surface normal with respect to the
corresponding zone. In this example, the normal is positive for
:math:`\alpha` and negative for :math:`\beta` (as indicated by the
arrows in the figure).

.. note::

  In this discussion, the term *interface* refers generically to an edge
  (2D) or face (3d) across which a flux may be calculated. The global
  ordering of the zones determines which zone actually creates the
  interface entity, thus determining its orientation.

FleCSI-SP captures this information in its connectivity graph with the
following conventions:

* *Zone-to-Interface* |br|
  The interface id is stored using a *one's complement* strategy
  depending on its orientation with respect to the zone, i.e., if the
  interface has a negative orientation, the id is stored as ~id (one's
  complement), otherwise, the raw id is stored. The Burton
  specialization provides the ``sign`` function to access interface ids:

  .. code-block:: cpp

    template<typename T>
    auto sign(T const & id) {
      static_assert(std::is_unsigned_v<T>, "id must use unsigned type");
      const bool neg{id>~id};
      return std::pair{1-2*neg, neg ? T(~id) : id};
    }

  The return values of this function consist of a *multiplier*
  (:math:`\pm 1`) and the raw id of the interface. These can be used to
  correctly capture the flux contribution across the interface:

  .. code-block:: cpp

    for(auto c: m.cells<owned>()) {
      for(auto e: m.edges(c)) {
        auto const [mul, id] = sign(e);
        q[c] += mul*f[id]; // quantity `q` is the sum of fluxes in `f`
      } // for
    } // for

* *Interface-to-Zone* |br|
  The zone ids at an interface are stored using the ordering convention
  that the zone for which the interface has a positive orientation is
  stored in the :math:`0^\mathrm{th}` place. This convention exploits
  the fact that an interface can only be connected to, at most, two
  zones.

.. vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 :
