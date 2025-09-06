function analyze_res(res1, res2, angle1, angle2)
    disp('--- 最终结果对比 ---');
    fprintf('文件3 (%d°) -> 厚度: %.2f μm, ν_p: %.1f, Γ: %.1f\n', angle1, res1.thk, res1.drudeParam(1), res1.drudeParam(2));
    fprintf('文件4 (%d°) -> 厚度: %.2f μm, ν_p: %.1f, Γ: %.1f\n', angle2, res2.thk, res2.drudeParam(1), res2.drudeParam(2));
    thk1 = res1.thk; thk2 = res2.thk;
    thk_diff_percent = abs(thk1-thk2)/mean([thk1,thk2])*100;
    fprintf('两个角度测得的厚度分别为 %.2f μm 和 %.2f μm，相对差异为 %.2f%%。\n', thk1, thk2, thk_diff_percent);
    if thk_diff_percent < 5
        fprintf('   -> 结论: 厚度一致性良好，结果可靠。\n');
        fprintf('   -> 推荐厚度值: %.2f ± %.2f μm\n', mean([thk1, thk2]), std([thk1, thk2]));
    else
        fprintf('   -> 结论: 厚度差异较大。\n');
    end
end
