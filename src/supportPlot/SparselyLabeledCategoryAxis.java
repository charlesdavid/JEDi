package supportPlot;

import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.axis.AxisState;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryTick;
import org.jfree.text.TextBlock;
import org.jfree.ui.RectangleEdge;

/**
 * A category axis that only writes the tick labels on a few of the ticks.
 */
public class SparselyLabeledCategoryAxis extends CategoryAxis
{
	private static final long serialVersionUID = 478725789943763302L;

	/**
	 * The number of ticks to label.
	 */
	private final int labeledTicks;

	/**
	 * Construct an axis without a label.
	 * 
	 * @param labeledTicks show only this many labeled ticks
	 */
	public SparselyLabeledCategoryAxis(int labeledTicks)
	{
		this.labeledTicks = labeledTicks;
	}

	/**
	 * Construct and axis with a label.
	 * 
	 * @param labeledTicks show only this many labeled ticks
	 * @param label        the axis label
	 */
	public SparselyLabeledCategoryAxis(int labeledTicks, String label)
	{
		super(label);
		this.labeledTicks = labeledTicks;
	}

	@Override
	@SuppressWarnings("unchecked")
	public List<CategoryTick> refreshTicks(Graphics2D g2, AxisState state, Rectangle2D dataArea, RectangleEdge edge)
	{
		List<CategoryTick> standardTicks = super.refreshTicks(g2, state, dataArea, edge);
		if (standardTicks.isEmpty())
			{
				return standardTicks;
			}
		int tickEvery = standardTicks.size() / labeledTicks;
		if (tickEvery < 1)
			{
				return standardTicks;
			}

		// Replace a few labels with blank ones
		List<CategoryTick> fixedTicks = new ArrayList<CategoryTick>(standardTicks.size());
		// Skip the first tick so your 45degree labels don't fall of the edge
		CategoryTick tick = standardTicks.get(0);
		fixedTicks.add(new CategoryTick(tick.getCategory(), new TextBlock(), tick.getLabelAnchor(), tick.getRotationAnchor(), tick.getAngle()));
		for (int i = 1; i < standardTicks.size(); i++)
			{
				tick = standardTicks.get(i);
				if (i % tickEvery == 0)
					{
						fixedTicks.add(tick);
					}
				else
					{
						fixedTicks.add(new CategoryTick(tick.getCategory(), new TextBlock(), tick.getLabelAnchor(), tick.getRotationAnchor(), tick.getAngle()));
					}
			}
		return fixedTicks;
	}
}
